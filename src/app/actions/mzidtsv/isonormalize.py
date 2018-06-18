import sys
from statistics import median, StatisticsError
from collections import OrderedDict

from app.dataformats import prottable as prottabledata
from app.readers import tsv as reader
from app.actions import isonormalizing


ISOQUANTRATIO_FEAT_ACC = '##isoquant_target_acc##'


def get_isobaric_ratios(psmfn, psmheader, channels, denom_channels, min_int,
                        targetfn, accessioncol, normalize, normratiofn):
    """Main function to calculate ratios for PSMs, peptides, proteins, genes.
    Can do simple ratios, median-of-ratios and median-centering
    normalization."""
    psm_or_feat_ratios = get_psmratios(psmfn, psmheader, channels,
                                       denom_channels, min_int, accessioncol)
    if normalize and normratiofn:
        normheader = reader.get_tsv_header(normratiofn)
        normratios = get_ratios_from_fn(normratiofn, normheader, channels)
        ch_medians = get_medians(channels, normratios, report=True)
        outratios = calculate_normalized_ratios(psm_or_feat_ratios, ch_medians,
                                                channels)
    elif normalize:
        flatratios = [[feat[ch] for ch in channels]
                      for feat in psm_or_feat_ratios]
        ch_medians = get_medians(channels, flatratios, report=True)
        outratios = calculate_normalized_ratios(psm_or_feat_ratios, ch_medians,
                                                channels)
    else:
        outratios = psm_or_feat_ratios
    # at this point, outratios look like:
    # [{ch1: 123, ch2: 456, ISOQUANTRATIO_FEAT_ACC: ENSG1244}, ]
    if accessioncol and targetfn:
        outratios = {x[ISOQUANTRATIO_FEAT_ACC]: x for x in outratios}
        return output_to_target_accession_table(targetfn, outratios, channels)
    elif not accessioncol and not targetfn:
        return paste_to_psmtable(psmfn, psmheader, outratios)
    elif accessioncol and not targetfn:
        # generate new table with accessions
        return ({(k if not k == ISOQUANTRATIO_FEAT_ACC
                  else prottabledata.HEADER_ACCESSION): v
                 for k, v in ratio.items()} for ratio in outratios)


def get_psmratios(psmfn, header, channels, denom_channels, min_int, acc_col):
    allfeats, feat_order, psmratios = {}, OrderedDict(), []
    for psm in reader.generate_tsv_psms(psmfn, header):
        ratios = calc_psm_ratios(psm, channels, denom_channels, min_int)
        # remove uninformative psms when adding to features
        if acc_col and (psm[acc_col] == '' or ';' in psm[acc_col] or
                        not {psm[q] for q in channels}.difference(
                            {'NA', None, False, ''})):
            continue
        elif acc_col:
            try:
                allfeats[psm[acc_col]].append(ratios)
            except KeyError:
                allfeats[psm[acc_col]] = [ratios]
            feat_order[psm[acc_col]] = 1
        else:
            psmquant = {ch: str(ratios[ix]) if ratios[ix] != 'NA' else 'NA'
                        for ix, ch in enumerate(channels)}
            psmquant[ISOQUANTRATIO_FEAT_ACC] = False
            psmratios.append(psmquant)
    if not acc_col:
        return psmratios
    else:
        outfeatures = []
        for feat in feat_order.keys():
            quants = allfeats[feat]
            outfeature = {ISOQUANTRATIO_FEAT_ACC: feat}
            outfeature.update(get_medians(channels, quants))
            outfeature.update(get_no_psms(channels, quants))
            outfeatures.append(outfeature)
    return outfeatures


def get_ratios_from_fn(fn, header, channels):
    ratios = []
    for feat in reader.generate_tsv_psms(fn, header):
        ratios.append([float(feat[ch]) if feat[ch] != 'NA' else 'NA'
                       for ch in channels])
    return ratios


def paste_to_psmtable(psmfn, header, ratios):
    # loop psms in psmtable, paste the outratios in memory
    for psm, ratio in zip(reader.generate_tsv_psms(psmfn, header), ratios):
        ratio.pop(ISOQUANTRATIO_FEAT_ACC)
        ratio = {'ratio_{}'.format(ch): val for ch, val in ratio.items()}
        psm.update(ratio)
        yield psm


def output_to_target_accession_table(targetfn, featratios, channels):
    #loop prottable, add ratios from dict, acc = key
    theader = reader.get_tsv_header(targetfn)
    acc_field = theader[0]
    for feat in reader.generate_tsv_proteins(targetfn, theader):
        try:
            quants = featratios[feat[acc_field]]
        except KeyError:
            quants = {ch: 'NA' for ch in channels}
            quants.update({get_no_psms_field(ch): 'NA' for ch in channels})
        else:
            quants.pop(ISOQUANTRATIO_FEAT_ACC)
        feat.update(quants)
        yield feat


def calc_psm_ratios(psm, channels, denom_channels, min_intensity):
    # set values below min_intensity to NA
    psm_intensity = {ch: float(psm[ch])
                     if psm[ch] != 'NA' and float(psm[ch]) > min_intensity
                     else 'NA' for ch in channels}
    denomvalues = [psm_intensity[ch] for ch in denom_channels
                   if psm_intensity[ch] != 'NA']
    if sum(denomvalues) == 0 or len(denomvalues) == 0:
        return ['NA'] * len(channels)
    denom = sum(denomvalues) / len(denomvalues)
    return [psm_intensity[ch] / denom
            if psm_intensity[ch] != 'NA' else 'NA' for ch in channels]


def get_medians(channels, ratios, report=False):
    ch_medians = {}
    for ix, channel in enumerate(channels):
        try:
            ch_medians[channel] = median([x[ix] for x in ratios
                                          if x[ix] != 'NA'])
        except StatisticsError:
            # channel is empty, common in protein quant but not in normalizing
            ch_medians[channel] = 'NA'
    if report:
        report = ('Channel intensity medians used for normalization:\n'
                  '{}'.format('\n'.join(['{} - {}'.format(ch, ch_medians[ch])
                                         for ch in channels])))
        sys.stdout.write(report)
    return ch_medians


def get_no_psms_field(quantfield):
    return '{}{}'.format(quantfield, prottabledata.HEADER_NO_PSMS_SUFFIX)


def get_no_psms(channels, ratios):
    ch_nopsms = {}
    for ix, channel in enumerate(channels):
        fieldname = get_no_psms_field(channel)
        ch_nopsms[fieldname] = len([x[ix] for x in ratios if x[ix] != 'NA'])
    return ch_nopsms


def calculate_normalized_ratios(ratios, ch_medians, channels):
    """Calculates ratios for PSM tables containing isobaric channels with
    raw intensities. Normalizes the ratios by median. NA values or values
    below min_intensity are excluded from the normalization."""
    outratios = []
    for quant in ratios:
        quant.update({ch: str(quant[ch] / ch_medians[ch])
                      if quant[ch] != 'NA' else 'NA' for ch in channels})
        outratios.append(quant)
    return outratios


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
