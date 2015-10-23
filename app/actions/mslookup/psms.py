from collections import OrderedDict

from app.readers import tsv as tsvreader
from app.readers import fasta as fastareader
DB_STORE_CHUNK = 100000


def create_psm_lookup(fn, fastafn, mapfn, header, pgdb, unroll=False,
                      specfncol=None):
    """Reads PSMs from file, stores them to a database backend in chunked PSMs.
    """
    store_proteins_descriptions(pgdb, fastafn, fn, mapfn, header)
    mzmlmap = pgdb.get_mzmlfile_map()
    sequences = {}
    for psm in tsvreader.generate_tsv_lines_multifile(fn, header):
        seq = tsvreader.get_psm_sequence(psm, unroll)
        sequences[seq] = 1
    pgdb.store_pepseqs(((seq,) for seq in sequences))
    pepseqmap = pgdb.get_peptide_seq_map()
    psms = []
    for rownr, psm in enumerate(tsvreader.generate_tsv_lines_multifile(fn, header)):
        specfn, psm_id, scan, seq, score = tsvreader.get_psm(psm, unroll, specfncol)
        if len(psms) % DB_STORE_CHUNK == 0:
            pgdb.store_psms(psms)
            psms = []
        psms.append({'rownr': rownr,
                     'psm_id': psm_id,
                     'seq': pepseqmap[seq],
                     'score': score,
                     'specfn': mzmlmap[specfn],
                     'scannr': scan,
                     })
    pgdb.store_psms(psms)
    pgdb.index_psms()
    store_psm_protein_relations(fn, header, pgdb)


def get_protein_gene_map(fastafn):
    gpmap = {}
    for protein, gene, symbol, desc in fastareader.get_proteins_genes(fastafn):
        gpmap[protein] = {'gene': gene, 'symbol': symbol, 'desc': desc}
    return gpmap


def store_proteins_descriptions(pgdb, fastafn, tsvfn, mapfn, header):
    if not fastafn: 
        proteins = {}
        for psm in tsvreader.generate_tsv_lines_multifile(tsvfn, header):
            proteins.update({x: 1 for x in 
                             tsvreader.get_proteins_from_psm(psm)}) 
        pgdb.store_proteins(((protein,) for protein in proteins[0].keys()))
    else:
        proteins, sequences, evidences = fastareader.get_proteins_for_db(
            fastafn)
        pgdb.store_proteins(proteins, evidences, sequences)
        if not mapfn:
            protein_descriptions = fastareader.get_proteins_descriptions(fastafn)
            pgdb.store_descriptions(protein_descriptions)
    if mapfn:
        gpmap = get_protein_gene_map(fastafn)
        store_gene_and_associated_id(gpmap)


def store_psm_protein_relations(fn, header, pgdb):
    """Reads PSMs from file, extracts their proteins and peptides and passes
    them to a database backend in chunks.
    """
    # TODO do we need an OrderedDict or is regular dict enough?
    # Sorting for psm_id useful?
    allpsms = OrderedDict()
    last_id, psmids_to_store = None, set()
    store_soon = False
    for psm in tsvreader.generate_tsv_lines_multifile(fn, header):
        psm_id, prots = tsvreader.get_pepproteins(psm)
        try:
            allpsms[psm_id].extend(prots)
        except KeyError:
            allpsms[psm_id] = prots
        if len(psmids_to_store) % DB_STORE_CHUNK == 0:
            store_soon = True
        if store_soon and last_id != psm_id:
            pgdb.store_peptides_proteins(allpsms, psmids_to_store)
            store_soon = False
            psmids_to_store = set()
        psmids_to_store.add(psm_id)
        last_id = psm_id
    pgdb.store_peptides_proteins(allpsms, psmids_to_store)
    pgdb.index_protein_peptides()
    return allpsms


