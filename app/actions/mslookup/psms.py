from app.readers import tsv as tsvreader
DB_STORE_CHUNK = 100000


def create_psm_lookup(fn, header, pgdb, unroll=False, specfncol=None):
    """Reads PSMs from file, stores them to a database backend in chunked PSMs.
    """
    mzmlmap = pgdb.get_mzmlfile_map()
    psms = []
    for rownr, psm in enumerate(tsvreader.generate_tsv_lines_multifile(fn, header)):
        specfn, psm_id, scan, seq, score = tsvreader.get_psm(psm, unroll, specfncol)
        if len(psms) % DB_STORE_CHUNK == 0:
            pgdb.store_psms(psms)
            psms = []
        psms.append({'rownr': rownr,
                     'psm_id': psm_id,
                     'seq': seq,
                     'score': score,
                     'specfn': mzmlmap[specfn],
                     'scannr': scan,
                     })
    pgdb.store_psms(psms)
    pgdb.index_psms()
