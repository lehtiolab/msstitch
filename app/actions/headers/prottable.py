from collections import OrderedDict

from app.actions.headers.base import (generate_general_header,
                                      generate_headerfields)
from app.dataformats import prottable as prottabledata


def generate_header(headerfields, oldheader=False):
    """Returns a header as a list, ready to write to TSV file"""
    fieldtypes = ['proteindata', 'probability', 'proteinfdr', 'proteinpep',
                  'precursorquant', 'isoquant']
    return generate_general_header(headerfields, fieldtypes,
                                   prottabledata.HEADER_PROTEIN, oldheader)


def get_prottable_headerfields(headertypes, lookup=False, poolnames=False):
    """Called by driver to generate headerfields object"""
    field_defs= {'isoquant': get_isoquant_fields(lookup, poolnames),
                 'precursorquant': get_precursorquant_fields(poolnames),
                 'probability': get_probability_fields(poolnames),
                 'proteindata': get_proteininfo_fields(poolnames),
                 'proteinfdr': get_proteinfdr_fields(poolnames),
                 'proteinpep': get_proteinpep_fields(poolnames),
                 }
    return generate_headerfields(headertypes, field_defs, poolnames)


def get_precursorquant_fields(poolnames=False):
    return {prottabledata.HEADER_AREA: poolnames}


def get_probability_fields(poolnames=False):
    return {prottabledata.HEADER_PROBABILITY: poolnames}


def get_proteinfdr_fields(poolnames=False):
    return {prottabledata.HEADER_QVAL: poolnames}


def get_proteinpep_fields(poolnames=False):
    return {prottabledata.HEADER_PEP: poolnames}


def get_proteininfo_fields(poolnames=False):
    """Returns header fields for protein (group) information.
    Some fields are shared between pools, others are specific
    for a pool"""
    allfields = OrderedDict()
    basefields = [prottabledata.HEADER_DESCRIPTION,
                  prottabledata.HEADER_COVERAGE,
                  prottabledata.HEADER_NO_PROTEIN,
                  ]
    poolfields = [prottabledata.HEADER_NO_UNIPEP,
                  prottabledata.HEADER_NO_PEPTIDE,
                  prottabledata.HEADER_NO_PSM,
                  ]
    for field in basefields:
        allfields[field] = False
    for field in poolfields:
        allfields[field] = poolnames
    return allfields


def get_isoquant_fields(pqdb=False, poolnames=False):
    """Returns a headerfield dict for isobaric quant channels. Channels are
    taken from DB and there isn't a pool-independent version of this yet"""
    quantheader = OrderedDict()
    for chan_name, amnt_psms_name in pqdb.get_isoquant_amountpsms_channels():
        quantheader[chan_name] = poolnames
        if amnt_psms_name:
            quantheader[amnt_psms_name] = poolnames
    return quantheader
