from collections import OrderedDict

from app.dataformats import mzidtsv as mzidtsvdata
from app.readers import tsv as tsvreader
from app.readers import fasta as fastareader
DB_STORE_CHUNK = 100000


def create_psm_lookup(fn, fastafn, header, pgdb, unroll=False,
                      specfncol=None, fastadelim=None, genefield=None):
    """Reads PSMs from file, stores them to a database backend in chunked PSMs.
    """
    proteins = store_proteins_descriptions(pgdb, fastafn, fn, header,
                                           fastadelim, genefield)
    mzmlmap = pgdb.get_mzmlfile_map()
    sequences = {}
    for psm in tsvreader.generate_tsv_psms(fn, header):
        seq = tsvreader.get_psm_sequence(psm, unroll)
        sequences[seq] = 1
    pgdb.store_pepseqs(((seq,) for seq in sequences))
    pepseqmap = pgdb.get_peptide_seq_map()
    psms = []
    for row, psm in enumerate(tsvreader.generate_tsv_psms(fn, header)):
        specfn, psm_id, specscanid, seq, score = tsvreader.get_psm(psm, unroll,
                                                             specfncol)
        if len(psms) % DB_STORE_CHUNK == 0:
            pgdb.store_psms(psms)
            psms = []
        psms.append({'rownr': row,
                     'psm_id': psm_id,
                     'seq': pepseqmap[seq],
                     'score': score,
                     'specfn': mzmlmap[specfn],
                     'spec_id': '{}_{}'.format(mzmlmap[specfn], specscanid),
                     })
    pgdb.store_psms(psms)
    pgdb.index_psms()
    store_psm_protein_relations(fn, header, pgdb, proteins)


def store_proteins_descriptions(pgdb, fastafn, tsvfn, header, fastadelim, genefield):
    if not fastafn:
        prots = {}
        for psm in tsvreader.generate_tsv_psms(tsvfn, header):
            prots.update({x: 1 for x in
                             tsvreader.get_proteins_from_psm(psm)})
        prots = [(protein,) for protein in prots.keys()]
        pgdb.store_proteins(prots)
    else:
        prots, seqs, desc, evids, ensgs, symbols = fastareader.get_proteins_for_db(
            fastafn, fastadelim, genefield)
        pgdb.store_fasta(prots, evids, seqs, desc, ensgs, symbols)
    return set([x[0] for x in prots])


def store_psm_protein_relations(fn, header, pgdb, proteins):
    """Reads PSMs from file, extracts their proteins and peptides and passes
    them to a database backend in chunks.
    """
    # TODO do we need an OrderedDict or is regular dict enough?
    # Sorting for psm_id useful?
    allpsms = OrderedDict()
    last_id, psmids_to_store = None, set()
    store_soon = False
    for psm in tsvreader.generate_tsv_psms(fn, header):
        psm_id, prots = tsvreader.get_pepproteins(psm)
        # TODO can this be removed permanently? 
        # Filter proteins to only include those that match the protein 
        # accessions in fasta so we get the correct names, filter out the badly annotated peptides
        # prots = [x for x in prots if x in proteins]
        try:
            # In case the PSMs are presented unrolled
            allpsms[psm_id].extend(prots)
        except KeyError:
            allpsms[psm_id] = prots
        if len(psmids_to_store) > DB_STORE_CHUNK:
            store_soon = True
        if store_soon and last_id != psm_id:
            pgdb.store_peptides_proteins(allpsms, psmids_to_store)
            store_soon = False
            psmids_to_store = set()
        psmids_to_store.add(psm_id)
        last_id = psm_id
    if len(psmids_to_store) > 0:
        pgdb.store_peptides_proteins(allpsms, psmids_to_store)
    pgdb.index_protein_peptides()
    return allpsms
