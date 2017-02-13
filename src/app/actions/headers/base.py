from collections import OrderedDict
from sqlite3 import OperationalError


def generate_general_header(headerfields, fieldtypes, firstfield,
                            oldheader=False):
    """From headerfield object, this generates a full header as a list, ready
    to write to a TSV file"""
    if not oldheader:
        header = [firstfield]
    else:
        header = oldheader[:]
    poolfields = OrderedDict()
    poolfields[None] = []  # Have non-pool/set columns come before pool-columns
    for fieldtype in fieldtypes:
        try:
            fields = headerfields[fieldtype]
        except KeyError:
            continue
        if type(fields) == list:
            header.extend(fields)
        else:
            for pool_field in fields.values():
                for pool, field in pool_field.items():
                    try:
                        poolfields[pool].append(field)
                    except KeyError:
                        poolfields[pool] = [field]
    if poolfields:
        for fields in poolfields.values():
            header.extend(fields)
    return header


def generate_headerfields(headertypes, allfield_defs, poolnames, db=False):
    """Returns a headerfields object (dict) which contains information on
    fields of the header, including optional pool names"""
    hfields = {}
    for fieldtype in headertypes:
        hfields[fieldtype] = OrderedDict()
        if fieldtype == 'isoquant':
            hfield_definitions = allfield_defs[fieldtype](db, poolnames)
        else:
            hfield_definitions = allfield_defs[fieldtype](poolnames)
        for fieldname, poolnames in hfield_definitions.items():
            hfields[fieldtype][fieldname] = get_header_field(fieldname,
                                                             poolnames)
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
    quantheader = OrderedDict()
    # FIXME when is a None database passed?
    if pqdb is None:
        return {}
    try:
        channels_psms = pqdb.get_isoquant_amountpsms_channels()
    except OperationalError:
        # FIXME what does this catch?
        return {}
    for chan_name, amnt_psms_name in channels_psms:
        quantheader[chan_name] = poolnames
        if amnt_psms_name:
            quantheader[amnt_psms_name] = poolnames
    return quantheader
