from collections import OrderedDict

from app.actions.headers.base import (generate_general_header,
                                      generate_headerfields)
from app.dataformats import prottable as prottabledata
from sqlite3 import OperationalError


def generate_header(headerfields, oldheader=False):
    """Returns a header as a list, ready to write to TSV file"""
    fieldtypes = ['proteindata', 'probability', 'proteinfdr', 'proteinpep',
                  'precursorquant', 'isoquant', 'bestpepscore']
    return generate_general_header(headerfields, fieldtypes,
                                   prottabledata.HEADER_PROTEIN, oldheader)


def get_prottable_headerfields(headertypes, lookup=False, poolnames=False):
    """Called by driver to generate headerfields object"""
    field_defs = {'isoquant': get_isoquant_fields,
                  'precursorquant': get_precursorquant_fields,
                  'probability': get_probability_fields,
                  'proteindata': get_proteininfo_fields,
                  'proteinfdr': get_proteinfdr_fields,
                  'proteinpep': get_proteinpep_fields,
                  'bestpepscore': get_bestpeptide_fields,
                  }
    return generate_headerfields(headertypes, field_defs, poolnames, lookup)


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
    basefields = [prottabledata.HEADER_GENE,
                  prottabledata.HEADER_ASSOCIATED,
                  prottabledata.HEADER_DESCRIPTION,
                  prottabledata.HEADER_COVERAGE,
                  prottabledata.HEADER_NO_PROTEIN,
                  prottabledata.HEADER_CONTENTPROT,
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
    if pqdb is None:
        return {}
    quantheader = OrderedDict()
    try:
        channels_psms = pqdb.get_isoquant_amountpsms_channels()
    except OperationalError:
        return {}
    for chan_name, amnt_psms_name in channels_psms:
        quantheader[chan_name] = poolnames
        if amnt_psms_name:
            quantheader[amnt_psms_name] = poolnames
    return quantheader


def get_bestpeptide_fields(poolnames=False):
    return {prottabledata.HEADER_QSCORE: poolnames}
