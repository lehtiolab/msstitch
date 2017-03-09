from app.actions.isonormalizing import get_medians
from app.dataformats import prottable as prottabledata
from app.actions.shared.pepprot_isoquant import base_add_isoquant_data


def get_isoquant_header(quantfields):
    baseheader = [prottabledata.HEADER_PROTEIN] + quantfields
    nopsms = [get_no_psms_field(qf) for qf in quantfields]
    return baseheader + nopsms


def get_no_psms_field(quantfield):
    return '{}{}'.format(quantfield, prottabledata.HEADER_NO_PSMS_SUFFIX)


def add_isoquant_data(proteins, quantproteins, quantacc, quantfields):
    """Runs through a protein table and adds quant data from ANOTHER protein
    table that contains that data."""
    for protein in base_add_isoquant_data(proteins, quantproteins,
                                          prottabledata.HEADER_PROTEIN,
                                          quantacc, quantfields):
        yield protein


def isobaric_quant_psms(psms, protcol, quantcols, normalize, normalize_table):
    """Runs through PSM table and uses medians of isobaric quant ratios of the
    PSMs per protein/peptide to create protein quantification output.
    Normalization can be applied either using the resulting ratios of the
    protein/peptide/genes or using a second list of features, passed as
    an iterable from the tsv reader.
    """
    features = get_isobaric_median_features(psms, protcol, quantcols)
    if normalize == 'median':
        if normalize_table:
            featratios = [[convert_to_float_or_na(f[x]) for x in quantcols]
                          for f in normalize_table]
        else:
            featratios = [[f[x] for x in quantcols] for f in features]
        ch_medians = get_medians(quantcols, featratios)
        for feat in features:
            feat.update({ch: str(feat[ch] / ch_medians[ch])
                         if feat[ch] != 'NA' else 'NA' for ch in quantcols})
    for feat in features:
        yield feat


def get_isobaric_median_features(psms, protcol, quantcols):
    protpsmquant, proteins = {}, []
    for psm in psms:
        if psm[protcol] == '' or ';' in psm[protcol]:
            continue
        elif not {psm[q] for q in quantcols}.difference({'NA',
                                                         None, False, ''}):
            continue
        try:
            protpsmquant[psm[protcol]].append([convert_to_float_or_na(psm[q])
                                               for q in quantcols])
        except KeyError:
            protpsmquant[psm[protcol]] = [[convert_to_float_or_na(psm[q])
                                           for q in quantcols]]
        proteins.append(psm[protcol])
    features = []
    for protein in proteins:
        quants = protpsmquant[protein]
        feature = {prottabledata.HEADER_PROTEIN: protein}
        feature.update(get_medians(quantcols, quants))
        feature.update(get_no_psms(quantcols, quants))
        features.append(feature)
    return features


def get_no_psms(channels, ratios):
    ch_nopsms = {}
    for ix, channel in enumerate(channels):
        fieldname = get_no_psms_field(channel)
        ch_nopsms[fieldname] = len([x[ix] for x in ratios if x[ix] != 'NA'])
    return ch_nopsms


def convert_to_float_or_na(value):
    try:
        return float(value)
    except ValueError:
        return 'NA'
