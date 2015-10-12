from app.actions.mzidtsv import confidencefilters as conffilt
from app.readers import tsv as tsvreader


def generate_psms(fn, header, confkey, conflvl, lower_is_better):
    for psm in tsvreader.generate_tsv_psms(fn, header):
        if conffilt.passes_filter(psm, conflvl, confkey, lower_is_better):
            yield psm
