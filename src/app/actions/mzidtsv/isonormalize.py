from app.readers import tsv as reader
from statistics import median


def calc_psm_ratios(psm, channels, denom_intensity):
    if denom_intensity == 0:
        return ['NA'] * len(channels)
    return [psm[ch] / denom_intensity
            if psm[ch] != 'NA' else 'NA' for ch in channels]


def calc_normalized(psm, channel, medians):
    return psm[channel] / medians[channel]


def get_normalized_ratios(psmfn, header, channels, denom_channels,
                          min_intensity):
    """Calculates ratios for PSM tables containing isobaric channels with
    raw intensities. Normalizes the ratios by median. NA values or values
    below min_intensity are excluded from the normalization."""
    ratios = []
    for psm in reader.generate_tsv_psms(psmfn, header):
        # set values below min_intensity to NA
        psmratios = {ch: float(psm[ch])
                     if psm[ch] != 'NA' and float(psm[ch]) > min_intensity
                     else 'NA' for ch in channels}
        denom = sum([psmratios[ch] for ch in denom_channels
                     if psmratios[ch] != 'NA']) / len(denom_channels)
        ratios.append(calc_psm_ratios(psmratios, channels, denom))
    ch_medians = get_medians(channels, ratios)
    for psm, psmratios in zip(reader.generate_tsv_psms(psmfn, header), ratios):
        psm.update({ch: str(psmratios[ix] / ch_medians[ch])
                    if psmratios[ix] != 'NA' else 'NA'
                    for ix, ch in enumerate(channels)})
        yield psm


def get_medians(channels, ratios):
    ch_medians = {}
    for ix, channel in enumerate(channels):
        ch_medians[channel] = median([x[ix] for x in ratios if x[ix] != 'NA'])
    return ch_medians
