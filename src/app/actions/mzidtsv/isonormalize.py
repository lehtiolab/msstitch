import sys

from app.actions.isonormalizing import get_medians
from app.dataformats import prottable as prottabledata
from app.readers import tsv as reader
from app.actions import isonormalizing


def get_isobaric_ratios(psmfn, psmheader, channels, denom_channels, min_int,
                        targetfn, accessioncol, normalize, normratiofn):
    # outputs:
    # PSM ratios
    # PSM ratios but normalized
    # protein median ratios
    # protein median normalized on itself
    # protein median normalized on another
    #
    psm_or_feat_ratios = get_psmratios(psmfn, psmheader, channels,
                                       denom_channels, min_int, accessioncol)
    if normalize and normratiofn:
        ch_medians = get_ch_medians(get_psmratios(normratiofn))
        outratios = generate_normalized_ratios(psm_or_feat_ratios, ch_medians)
    elif normalize:
        psm_or_feat_ratios = [x for x in psm_or_feat_ratios]
        ch_medians = get_ch_medians(channels, psmratios)
        outratios = generate_normalized_ratios(psm_or_feat_ratios, ch_medians)
    else:
        outratios = [x for x in psm_or_feat_ratios]
    if accessioncol and targetfn:
        output_to_target_accession_table(outratios, targetfn)
    elif targetfn == psmfn:
        return paste_to_psmtable(psmfn, psmheader, outratios)
    else:
        return outratios


def get_psmratios(psmfn, header, channels, denom_channels, min_int, acc_col):
    for psm in reader.generate_tsv_psms(psmfn, header):
        ratios = calc_psm_ratios(psm, channels, denom_channels, min_int)
        yield {ch: str(ratios[ix]) if ratios[ix] != 'NA' else 'NA'
               for ix, ch in enumerate(channels)}


def paste_to_psmtable(psmfn, header, ratios):
    # loop psms in psmtable, paste the outratios in memory
    for psm, ratio in zip(reader.generate_tsv_psms(psmfn, header), ratios):
        psm.update(ratio)
        yield psm


def output_to_target_accession_table():
    #loop prottable, add ratios from dict, acc = key
    pass


def calc_psm_ratios(psm, channels, denom_channels, min_intensity):
    # set values below min_intensity to NA
    psm_intensity = {ch: float(psm[ch])
                     if psm[ch] != 'NA' and float(psm[ch]) > min_intensity
                     else 'NA' for ch in channels}
    denomvalues = [psm_intensity[ch] for ch in denom_channels
                   if psm_intensity[ch] != 'NA']
    denom = sum(denomvalues) / len(denomvalues)
    if denom == 0:
        return ['NA'] * len(channels)
    return [psm_intensity[ch] / denom
            if psm_intensity[ch] != 'NA' else 'NA' for ch in channels]


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
    ch_medians = isonormalizing.get_medians(channels, ratios)
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
