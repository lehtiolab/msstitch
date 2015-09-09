from collections import OrderedDict

from app.actions.headers.base import (generate_general_header,
                                      generate_headerfields)
from app.dataformats import peptable as peptabledata


def get_pepquant_header(oldheader, isobqfieldmap=False, precurqfield=False):
    header = oldheader[:]
    if isobqfieldmap:
        for field, medianfield in isobqfieldmap.items():
            header = [medianfield if x == field else x for x in header]
    if precurqfield:
        header = [peptabledata.HEADER_AREA if x == precurqfield
                  else x for x in header]
    peptable_header = [peptabledata.HEADER_LINKED_PSMS]
    ix = header.index(peptabledata.HEADER_PEPTIDE)
    return header[:ix] + peptable_header + header[ix:]


def generate_header(headerfields, oldheader=False):
    """Returns a header as a list, ready to write to TSV file"""
    fieldtypes = ['peptidefdr', 'peptidepep',
                  'precursorquant', 'isoquant']
    return generate_general_header(headerfields, fieldtypes,
                                   peptabledata.HEADER_PEPTIDE, oldheader)


def get_prottable_headerfields(headertypes, lookup=False, poolnames=False):
    """Called by driver to generate headerfields object"""
    field_defs = {'isoquant': get_isoquant_fields(lookup, poolnames),
                  'precursorquant': get_precursorquant_fields(poolnames),
                  'peptidefdr': get_peptidefdr_fields(poolnames),
                  'peptidepep': get_peptidepep_fields(poolnames),
                  }
    return generate_headerfields(headertypes, field_defs, poolnames)


def get_precursorquant_fields(poolnames=False):
    return {peptabledata.HEADER_AREA: poolnames}


def get_peptidefdr_fields(poolnames=False):
    return {peptabledata.HEADER_QVAL: poolnames}


def get_peptidepep_fields(poolnames=False):
    return {peptabledata.HEADER_PEP: poolnames}


#def get_proteininfo_fields(poolnames=False):
#    """Returns header fields for protein (group) information.
#    Some fields are shared between pools, others are specific
#    for a pool"""
#    allfields = OrderedDict()
#    basefields = [prottabledata.HEADER_DESCRIPTION,
#                  prottabledata.HEADER_COVERAGE,
#                  prottabledata.HEADER_NO_PROTEIN,
#                  ]
#    poolfields = [prottabledata.HEADER_NO_UNIPEP,
#                  prottabledata.HEADER_NO_PEPTIDE,
#                  prottabledata.HEADER_NO_PSM,
#                  ]
#    for field in basefields:
#        allfields[field] = False
#    for field in poolfields:
#        allfields[field] = poolnames
#    return allfields


def get_isoquant_fields(pqdb=False, poolnames=False):
    """Returns a headerfield dict for isobaric quant channels. Channels are
    taken from DB and there isn't a pool-independent version of this yet"""
    quantheader = OrderedDict()
    for chan_name in pqdb.get_isoquant_channels():
        quantheader[chan_name] = poolnames
    return quantheader
