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


def get_no_psms(channels, ratios):
    ch_nopsms = {}
    for ix, channel in enumerate(channels):
        fieldname = get_no_psms_field(channel)
        ch_nopsms[fieldname] = len([x[ix] for x in ratios if x[ix] != 'NA'])
    return ch_nopsms
