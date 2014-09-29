def passes_filter(psm, threshold, confkey, lower_is_better):
    if psm[confkey] in ['NA', '', None, False]:
        return False
    lower = float(psm[confkey]) < float(threshold)
    return lower == lower_is_better
