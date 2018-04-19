from numpy import transpose
from collections import OrderedDict
from sqlite3 import OperationalError

from app.dataformats.prottable import HEADER_NO_PSMS_SUFFIX

def generate_general_header(headerfields, fieldtypes, firstfield, oldheader,
                            group_by_field):
    """From headerfield object, this generates a full header as a list, ready
    to write to a TSV file
    E.g:
    headerfield = {precusroquant: {HEADER_AREA: OD([(set1, set1_HEAD), (set2,
    set2_HEAD), etc])}}"""
    if not oldheader:
        header = [firstfield]
    else:
        header = oldheader[:]
    poolfields = OrderedDict()
    poolfields[None] = []  # Have non-pool/set columns come before pool-columns
    if group_by_field:
        header.extend(poolfields[None])
    for fieldtype in fieldtypes:
        try:
            fields = headerfields[fieldtype]
        except KeyError:
            continue
        if type(fields) == list:
            header.extend(fields)
        elif group_by_field:
            pfmatrix = [list(x.values()) for k, x in fields.items()
                        if not HEADER_NO_PSMS_SUFFIX in k]
            header.extend([x for y in transpose(pfmatrix) for x in y])
            if fieldtype == 'isoquant':
                pfmatrix = [list(x.values()) for k, x in fields.items()
                            if HEADER_NO_PSMS_SUFFIX in k]
                header.extend([x for y in transpose(pfmatrix) for x in y])
        else:
            for pool_field in fields.values():
                for pool, field in pool_field.items():
                    try:
                        poolfields[pool].append(field)
                    except KeyError:
                        poolfields[pool] = [field]
    if poolfields and not group_by_field:
        for fields in poolfields.values():
            header.extend(fields)
    return header


def generate_headerfields(headertypes, allfield_defs, poolnames, db=False):
    """Returns a headerfields object (dict) which contains information on
    fields of the header, including optional pool names"""
    hfields = {}
    for fieldtype in headertypes:
        if fieldtype == 'isoquant':
            continue
        hfield_definitions = allfield_defs[fieldtype](poolnames)
        hfields[fieldtype] = OrderedDict()
        for fieldname, poolnames in hfield_definitions.items():
            hfields[fieldtype][fieldname] = get_header_field(fieldname,
                                                             poolnames)
    if 'isoquant' in headertypes:
        hfield_definitions = allfield_defs['isoquant'](db, poolnames)
        hfields['isoquant'] = OrderedDict()
        for poolname in poolnames:
            for fieldname in hfield_definitions:
                hfields['isoquant'][fieldname] = get_header_field(
                    fieldname, poolnames)
    return hfields


def get_header_field(field, poolnames=False):
    if poolnames:
        return OrderedDict([(pool, '{}_{}'.format(pool, field))
                            for pool in poolnames])
    else:
        return {None: field}


def get_isoquant_fields(pqdb=False, poolnames=False):
    """Returns a headerfield dict for isobaric quant channels. Channels are
    taken from DB and there isn't a pool-independent version of this yet"""
    # FIXME when is a None database passed?
    if pqdb is None:
        return {}
    try:
        channels_psms = pqdb.get_isoquant_amountpsms_channels()
    except OperationalError:
        # FIXME what does this catch?
        return {}
    quantheader, psmsheader = OrderedDict(), OrderedDict()
    for chan_name, amnt_psms_name in channels_psms:
        quantheader[chan_name] = poolnames
        if amnt_psms_name:
            psmsheader[amnt_psms_name] = poolnames
    quantheader.update(psmsheader)
    return quantheader
