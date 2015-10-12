from app.readers import tsv as tsvreader
from app.readers import fasta as fastareader
DB_STORE_CHUNK = 100000


def create_psm_lookup(fn, fastafn, header, pgdb, unroll=False, specfncol=None, 
                      proteinfield=False):
    """Reads PSMs from file, stores them to a database backend in chunked PSMs.
    """
    store_proteins_descriptions(pgdb, fastafn, fn, header, proteinfield)
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


def store_proteins_descriptions(pgdb, fastafn, tsvfn, header,
                                proteinfield=False):
    if fastafn:
        proteins, sequences, evidences = fastareader.get_proteins_for_db(
            fastafn)
        pgdb.store_proteins(proteins, evidences, sequences)
        protein_descriptions = fastareader.get_proteins_descriptions(fastafn)
        pgdb.store_descriptions(protein_descriptions)
    else:
        proteins = {}
        for psm in tsvreader.generate_tsv_lines_multifile(tsvfn, header):
            proteins.update({x: 1 for x in
                             tsvreader.get_proteins_from_psm(psm,
                                                             proteinfield)})
        pgdb.store_proteins(((protein,) for protein in proteins.keys()))


