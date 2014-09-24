def passes_filter(psm, threshold, confkey, lower_is_better):
    lower = float(psm[confkey]) < float(threshold)
    return lower == lower_is_better
