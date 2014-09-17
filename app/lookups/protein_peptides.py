MZIDTSV_PEP_COL = 9
MZIDTSV_PROT_COL = 10
DB_STORE_CHUNK = 10000

from app import sqlite
from app.readers import tsv as tsvreader


def create_protein_pep_lookup(fn, unroll=False):
    """Reads PSMs from file, extracts their proteins and peptides and passes
    them to a database backend in chunked PSMs.
    """
    ppdb = sqlite.ProteinPeptideDB()
    ppdb.create_ppdb()
    psms = tsvreader.get_mzidtsv_lines_scannr_specfn(fn)
    count, last_id, peptides_proteins = 0, None, {}
    for line, (scan, specfn) in psms:
        count += 1
        pep_id, seq, prots = tsvreader.get_peptide_proteins(line, specfn,
                                                            scan, unroll)
        if count >= DB_STORE_CHUNK and pep_id != last_id:
            ppdb.store_peptides_proteins(peptides_proteins)
            peptides_proteins = {}
        try:
            peptides_proteins[pep_id]['proteins'].append(prots)
        except KeyError:
            peptides_proteins[pep_id] = {'scan_nr': scan, 'specfn': specfn,
                                         'seq': seq, 'proteins': prots}
        last_id = pep_id
    ppdb.store_peptides_proteins(peptides_proteins)
    ppdb.index()
    return ppdb.fn
