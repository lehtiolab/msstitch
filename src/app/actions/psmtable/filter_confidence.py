from app.readers import tsv as tsvreader


def filter_psms(psms, confkey, conflvl, lower_is_better):
    for psm in psms:
        if passes_filter(psm, conflvl, confkey, lower_is_better):
            yield psm


def passes_filter(psm, threshold, confkey, lower_is_better):
    try:
        confval = float(psm[confkey])
    except (TypeError, ValueError):
        return False
    else:
        lower = confval < float(threshold)
        return lower == lower_is_better
