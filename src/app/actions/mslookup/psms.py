from collections import OrderedDict

from app.dataformats import mzidtsv as mzidtsvdata
from app.readers import tsv as tsvreader
from app.readers import fasta as fastareader
DB_STORE_CHUNK = 100000


def create_psm_lookup(fn, fastafn, mapfn, header, pgdb, unroll=False,
                      specfncol=None, decoy=False,
                      fastadelim=None, genefield=None):
    """Reads PSMs from file, stores them to a database backend in chunked PSMs.
    """
    proteins = store_proteins_descriptions(pgdb, fastafn, fn, mapfn, header,
                                           decoy, fastadelim, genefield)
    mzmlmap = pgdb.get_mzmlfile_map()
    sequences = {}
    for psm in tsvreader.generate_tsv_psms(fn, header):
        seq = tsvreader.get_psm_sequence(psm, unroll)
        sequences[seq] = 1
    pgdb.store_pepseqs(((seq,) for seq in sequences))
    pepseqmap = pgdb.get_peptide_seq_map()
    psms = []
    for row, psm in enumerate(tsvreader.generate_tsv_psms(fn, header)):
        specfn, psm_id, scan, seq, score = tsvreader.get_psm(psm, unroll,
                                                             specfncol)
        if len(psms) % DB_STORE_CHUNK == 0:
            pgdb.store_psms(psms)
            psms = []
        psms.append({'rownr': row,
                     'psm_id': psm_id,
                     'seq': pepseqmap[seq],
                     'score': score,
                     'specfn': mzmlmap[specfn],
                     'scannr': scan,
                     'spec_id': '{}_{}'.format(mzmlmap[specfn], scan),
                     })
    pgdb.store_psms(psms)
    pgdb.index_psms()
    store_psm_protein_relations(fn, header, pgdb, proteins)


def get_protein_gene_map(mapfn, proteins, decoy):
    gpmap = {}
    for protein, gene, symbol, desc in fastareader.get_proteins_genes(mapfn):
        if protein in proteins:
            protein = proteins[protein][0]
            if decoy:
                symbol = '{}{}'.format(mzidtsvdata.DECOY_PREFIX, symbol)
                gene = '{}{}'.format(mzidtsvdata.DECOY_PREFIX, gene)
                desc = None
            gpmap[protein] = {'gene': gene, 'symbol': symbol, 'desc': desc}
    return gpmap


def store_proteins_descriptions(pgdb, fastafn, tsvfn, mapfn, header, decoy,
                                fastadelim, genefield):
    if not fastafn:
        proteins = {}
        for psm in tsvreader.generate_tsv_psms(tsvfn, header):
            proteins.update({x: 1 for x in
                             tsvreader.get_proteins_from_psm(psm)})
        proteins = [(protein,) for protein in proteins.keys()]
        pgdb.store_proteins(proteins)
    else:
        proteins, sequences, evidences = fastareader.get_proteins_for_db(
            fastafn)
        proteins = [x for x in proteins]
        pgdb.store_proteins(proteins, evidences, sequences)
        if not mapfn:
            associations = fastareader.get_proteins_genes(fastafn, fastadelim,
                                                          genefield)
            genes, descriptions = [], []
            for assoc in associations:
                genes.append((assoc[1], assoc[0]))
                descriptions.append((assoc[0], assoc[3]))
            pgdb.store_descriptions(descriptions)
            pgdb.store_genes(genes)
    if mapfn:
        proteins_with_versions = {}
        for protein in proteins:
            proteins_with_versions[protein[0].split('.')[0]] = protein
        if decoy:
            mod = fastareader.get_decoy_mod_string(proteins[0][0])
            proteins_with_versions = {k.replace(mod, ''): v for k, v
                                      in proteins_with_versions.items()}
        gpmap = get_protein_gene_map(mapfn, proteins_with_versions, decoy)
        pgdb.store_gene_and_associated_id(gpmap)
    return set([x[0] for x in proteins])


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
        prots = [x for x in prots if x in proteins]
        try:
            # In case the PSMs are presented unrolled
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
    if len(psmids_to_store) > 0:
        pgdb.store_peptides_proteins(allpsms, psmids_to_store)
    pgdb.index_protein_peptides()
    return allpsms
