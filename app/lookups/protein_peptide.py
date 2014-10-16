MZIDTSV_PEP_COL = 9
MZIDTSV_PROT_COL = 10
DB_STORE_CHUNK = 500000

from app.sqlite import ProteinGroupDB
from app.readers import tsv as tsvreader
from app.preparation.mzidtsv import confidencefilters as conffilt


def create_protein_pep_lookup(fn, header, confkey, conflvl, lower_is_better,
                              unroll=False, fastafn=None, evidence_lvl=False):
    """Reads PSMs from file, extracts their proteins and peptides and passes
    them to a database backend in chunked PSMs.
    """
    pgdb = ProteinGroupDB()
    pgdb.create_pgdb()
    rownr, last_id, peptides_proteins = 0, None, {}
    store_soon = False
    for psm in tsvreader.generate_tsv_psms(fn, header):
        if not conffilt.passes_filter(psm, conflvl, confkey, lower_is_better):
            rownr += 1
            continue
        specfn, scan, seq, score, prots = tsvreader.get_pepproteins(
            psm, unroll)
        psm_id = '{0}_{1}'.format(specfn, scan)
        if peptides_proteins and len(peptides_proteins) % DB_STORE_CHUNK == 0:
            store_soon = True
        if store_soon and last_id != psm_id:
            pgdb.store_peptides_proteins(peptides_proteins)
            store_soon = False
            peptides_proteins = {}
        try:
            peptides_proteins[psm_id]['proteins'].extend(prots)
            peptides_proteins[psm_id]['rows'].append(rownr)
        except KeyError:
            peptides_proteins[psm_id] = {
                                        'rows': [rownr],                                       
                                        'seq': seq, 
                                        'proteins': prots,
                                        'score': score,
                                        }
        last_id = psm_id
        rownr += 1
    pgdb.store_peptides_proteins(peptides_proteins)
    pgdb.index_protein_peptides()
    return pgdb
