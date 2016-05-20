from app.actions.mzidtsv.isonormalize import get_medians
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


def isobaric_quant_psms(psms, protcol, quantcols):
    """Runs through PSM table and uses medians of isobaric quant values of the
    PSMs per protein to create protein quantification output. Combine with
    median normalization for more accurate results"""
    protpsmquant = {}
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
    for protein, quants in protpsmquant.items():
        outprotein = {prottabledata.HEADER_PROTEIN: protein}
        outprotein.update(get_medians(quantcols, quants))
        outprotein.update(get_no_psms(quantcols, quants))
        yield outprotein


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
