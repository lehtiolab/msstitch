from app.dataformats import mzidtsv as header


def filter_psms_conf(psms, confkey, conflvl, lower_is_better):
    threshold = float(conflvl)
    for psm in psms:
        try:
            confval = float(psm[confkey])
        except (TypeError, ValueError):
            pass
        else:
            lower = confval < threshold
            if lower == lower_is_better:
                yield psm


def filter_psms_remove_set(psms, setnames):
    biosets = set(setnames)
    for psm in psms:
        if psm[header.HEADER_SETNAME] not in biosets:
            yield psm
