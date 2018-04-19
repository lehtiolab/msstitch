from collections import OrderedDict

from app.readers import tsv
from app.actions.headers.base import (generate_general_header,
                                      generate_headerfields,
                                      get_isoquant_fields)
from app.dataformats import peptable as peptabledata
from app.dataformats import mzidtsv as mzidtsvdata


def switch_psm_to_peptable_fields(oldheader):
    """Returns a dict map with old to new header fields"""
    return {old: new for old, new in zip([mzidtsvdata.HEADER_PEPTIDE,
                                          mzidtsvdata.HEADER_PROTEIN,
                                          mzidtsvdata.HEADER_PEPTIDE_Q,
                                          mzidtsvdata.HEADER_PEPTIDE_PEP],
                                         [peptabledata.HEADER_PEPTIDE,
                                          peptabledata.HEADER_PROTEINS,
                                          peptabledata.HEADER_QVAL,
                                          peptabledata.HEADER_PEP])}


def get_psm2pep_header(oldheader, isobq_pattern=False, precurqfield=False):
    header = oldheader[:]
    if isobq_pattern:
        isocols = tsv.get_columns_by_pattern(header, isobq_pattern)
        for col in isocols:
            header.pop(header.index(col))
    if precurqfield:
        header = [peptabledata.HEADER_AREA if x == precurqfield
                  else x for x in header]
    peptable_header = [peptabledata.HEADER_LINKED_PSMS]
    ix = header.index(mzidtsvdata.HEADER_PEPTIDE)
    header = header[:ix] + peptable_header + header[ix:]
    switch_map = switch_psm_to_peptable_fields(header)
    return [switch_map[field] if field in switch_map else field
            for field in header]


def get_linear_model_header(oldheader):
    header = oldheader[:]
    try:
        ix = header.index(peptabledata.HEADER_PEP) + 1
    except ValueError:
        ix = header.index(peptabledata.HEADER_QVAL) + 1
    return header[:ix] + [peptabledata.HEADER_QVAL_MODELED] + header[ix:]


def generate_header(headerfields, oldheader, group_by_field):
    """Returns a header as a list, ready to write to TSV file"""
    fieldtypes = ['peptidefdr', 'peptidepep', 'nopsms', 'proteindata',
                  'precursorquant', 'isoquant']
    return generate_general_header(headerfields, fieldtypes,
                                   peptabledata.HEADER_PEPTIDE, oldheader,
                                   group_by_field)


def get_peptable_headerfields(headertypes, lookup=False, poolnames=False):
    """Called by driver to generate headerfields object"""
    field_defs = {'isoquant': get_isoquant_fields,
                  'precursorquant': get_precursorquant_fields,
                  'peptidefdr': get_peptidefdr_fields,
                  'peptidepep': get_peptidepep_fields,
                  'proteindata': get_proteininfo_fields,
                  }
    return generate_headerfields(headertypes, field_defs, poolnames, lookup)


def get_precursorquant_fields(poolnames=False):
    return {peptabledata.HEADER_AREA: poolnames}


def get_peptidefdr_fields(poolnames=False):
    return {peptabledata.HEADER_QVAL: poolnames}


def get_peptidepep_fields(poolnames=False):
    return {peptabledata.HEADER_PEP: poolnames}


def get_proteininfo_fields(poolnames=False):
    """Returns header fields for protein (group) information."""
    allfields = OrderedDict()
    basefields = [peptabledata.HEADER_PROTEINS,
                  peptabledata.HEADER_GENES,
                  peptabledata.HEADER_ASSOCIATED,
                  peptabledata.HEADER_DESCRIPTIONS,
                  peptabledata.HEADER_COVERAGES,
                  peptabledata.HEADER_NO_CONTENTPROTEINS,
                  ]
    for field in basefields:
        allfields[field] = False
    allfields[peptabledata.HEADER_NO_PSM] = poolnames
    return allfields
