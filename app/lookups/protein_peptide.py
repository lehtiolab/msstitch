MZIDTSV_PEP_COL = 9
MZIDTSV_PROT_COL = 10
DB_STORE_CHUNK = 500000

from app.sqlite import ProteinGroupDB
from app.readers import tsv as tsvreader
from app.preparation.mzidtsv import confidencefilters as conffilt


def create_protein_pep_lookup(fn, header, confkey, conflvl, lower_is_better,
                              unroll=False):
    """Reads PSMs from file, extracts their proteins and peptides and passes
    them to a database backend in chunked PSMs.
    """
    pgdb = ProteinGroupDB()
    pgdb.create_pgdb()
    rownr, last_rownr, peptides_proteins = 0, None, {}
    for psm in tsvreader.generate_tsv_psms(fn, header):
        if not conffilt.passes_filter(psm, conflvl, confkey, lower_is_better):
            rownr += 1
            continue
        specfn, scan, seq, score, prots = tsvreader.get_pepproteins(
            psm, unroll)
        if rownr % DB_STORE_CHUNK == 0 and rownr != last_rownr:
            pgdb.store_peptides_proteins(peptides_proteins)
            peptides_proteins = {}
        try:
            peptides_proteins[rownr]['proteins'].extend(prots)
        except KeyError:
            peptides_proteins[rownr] = {'seq': seq, 
                                        'proteins': prots,
                                        'score': score,
                                        }
        last_rownr = rownr
        rownr += 1
    pgdb.store_peptides_proteins(peptides_proteins)
    pgdb.index_protein_peptides()
    return pgdb
