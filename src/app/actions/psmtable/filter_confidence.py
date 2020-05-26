from app.readers import tsv as tsvreader


def filter_psms(psms, confkey, conflvl, lower_is_better):
    for psm in psms:
        if passes_filter(psm, conflvl, confkey, lower_is_better):
            yield psm


def passes_filter(psm, threshold, confkey, lower_is_better):
    if psm[confkey] in ['NA', '', None, False]:
        return False
    lower = float(psm[confkey]) < float(threshold)
    return lower == lower_is_better
