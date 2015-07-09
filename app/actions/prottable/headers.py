from collections import OrderedDict

from app.dataformats import prottable as prottabledata


def get_headerfields(headertypes, lookup=False, poolnames=False):
    hfields = {}
    type_functions = {'isoquant': get_isoquant_fields,
                      'precursorquant': get_precursorquant_fields,
                      'probability': get_probability_fields,
                      'proteindata': get_proteininfo_fields,
                      'proteinfdr': get_proteinfdr_fields,
                      'proteinpep': get_proteinpep_fields,
                      }
    for fieldtype in headertypes:
        hfields[fieldtype] = type_functions[fieldtype](lookup, poolnames)
    return hfields


def get_header_field(field, poolnames=False):
    if poolnames:
        return OrderedDict([(pool, '{}_{}'.format(pool, field))
                            for pool in poolnames])
    else:
        return {None: field}


def get_precursorquant_fields(pqdb=False, poolnames=False):
    return {prottabledata.HEADER_AREA:
            get_header_field(prottabledata.HEADER_AREA, poolnames)}


def get_probability_fields(pqdb=False, poolnames=False):
    return {prottabledata.HEADER_PROBABILITY:
            get_header_field(prottabledata.HEADER_PROBABILITY, poolnames)}


def get_proteinfdr_fields(pqbd=False, poolnames=False):
    return {prottabledata.HEADER_QVAL:
            get_header_field(prottabledata.HEADER_QVAL, poolnames)}


def get_proteinpep_fields(pqbd=False, poolnames=False):
    return {prottabledata.HEADER_PEP:
            get_header_field(prottabledata.HEADER_PEP, poolnames)}


def get_proteininfo_fields(pqdb=False, poolnames=False):
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
        allfields[field] = get_header_field(field)
    for field in poolfields:
        allfields[field] = get_header_field(field, poolnames)
    return allfields


def get_isoquant_fields(pqdb=False, poolnames=False):
    """Returns a headerfield dict for isobaric quant channels. Channels are
    taken from DB and there isn't a pool-independent version of this yet"""
    quantheader = OrderedDict()
    for chan_name, amnt_psms_name in pqdb.get_isoquant_amountpsms_channels():
        quantheader[chan_name] = get_header_field(chan_name, poolnames)
        quantheader[amnt_psms_name] = get_header_field(amnt_psms_name,
                                                       poolnames)
    return quantheader


def generate_header(headerfields, oldheader=False):
    if not oldheader:
        header = [prottabledata.HEADER_PROTEIN]
    else:
        header = oldheader[:]
    for fieldtype in ['proteindata', 'probability', 'precursorquant',
                      'isoquant']:
        try:
            fields = headerfields[fieldtype]
        except KeyError:
            continue
        if type(fields) == list:
            header.extend(fields)
        else:
            for pools in fields.values():
                header.extend(pools.values())
    return header
