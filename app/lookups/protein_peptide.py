MZIDTSV_PEP_COL = 9
MZIDTSV_PROT_COL = 10
DB_STORE_CHUNK = 10000

from app import sqlite
from app.readers import tsv as tsvreader
from app.preparation.mzidtsv import confidencefilters as conffilt


def create_protein_pep_lookup(fn, header, confkey, conflvl, lower_is_better,
                              unroll=False):
    """Reads PSMs from file, extracts their proteins and peptides and passes
    them to a database backend in chunked PSMs.
    """
    ppdb = sqlite.ProteinPeptideDB()
    ppdb.create_ppdb()
    count, last_id, peptides_proteins = 0, None, {}
    for psm in tsvreader.generate_tsv_psms(fn, header):
        if not conffilt.passes_filter(psm, conflvl, confkey, lower_is_better):
            continue
        count += 1
        specfn, scan, psm_id, seq, prots = tsvreader.get_pepproteins(
            psm, unroll)
        if count >= DB_STORE_CHUNK and psm_id != last_id:
            ppdb.store_peptides_proteins(peptides_proteins)
            peptides_proteins = {}
        try:
            peptides_proteins[psm_id]['proteins'].extend(prots)
        except KeyError:
            peptides_proteins[psm_id] = {'scan_nr': scan, 'specfn': specfn,
                                         'seq': seq, 'proteins': prots}
        last_id = psm_id
    ppdb.store_peptides_proteins(peptides_proteins)
    ppdb.index()
    return ppdb.fn
