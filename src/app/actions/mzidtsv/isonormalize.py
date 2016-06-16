import sys
from statistics import median, StatisticsError

from app.readers import tsv as reader


def calc_psm_ratios(psm, channels, denom_channels, min_intensity):
    # set values below min_intensity to NA
    psm_intensity = {ch: float(psm[ch])
                     if psm[ch] != 'NA' and float(psm[ch]) > min_intensity
                     else 'NA' for ch in channels}
    denom = sum([psm_intensity[ch] for ch in denom_channels
                 if psm_intensity[ch] != 'NA']) / len(denom_channels)
    if denom == 0:
        return ['NA'] * len(channels)
    return [psm_intensity[ch] / denom
            if psm_intensity[ch] != 'NA' else 'NA' for ch in channels]


def calc_normalized(psm, channel, medians):
    return psm[channel] / medians[channel]


def get_normalized_ratios(psmfn, header, channels, denom_channels,
                          min_intensity, second_psmfn, secondheader):
    """Calculates ratios for PSM tables containing isobaric channels with
    raw intensities. Normalizes the ratios by median. NA values or values
    below min_intensity are excluded from the normalization."""
    ratios = []
    if second_psmfn is not None:
        median_psmfn = second_psmfn
        medianheader = secondheader
    else:
        median_psmfn = psmfn
        medianheader = header
    for psm in reader.generate_tsv_psms(median_psmfn, medianheader):
        ratios.append(calc_psm_ratios(psm, channels, denom_channels,
                                      min_intensity))
    ch_medians = get_medians(channels, ratios)
    report = ('Channel intensity medians used for normalization:\n'
              '{}'.format('\n'.join(['{} - {}'.format(ch, ch_medians[ch])
                                     for ch in channels])))
    sys.stdout.write(report)
    for psm in reader.generate_tsv_psms(psmfn, header):
        psmratios = calc_psm_ratios(psm, channels, denom_channels,
                                    min_intensity)
        psm.update({ch: str(psmratios[ix] / ch_medians[ch])
                    if psmratios[ix] != 'NA' else 'NA'
                    for ix, ch in enumerate(channels)})
        yield psm


def get_medians(channels, ratios):
    ch_medians = {}
    for ix, channel in enumerate(channels):
        try:
            ch_medians[channel] = median([x[ix] for x in ratios
                                          if x[ix] != 'NA'])
        except StatisticsError:
            # channel is empty, common in protein quant but not in normalizing
            ch_medians[channel] = 'NA'
    return ch_medians
