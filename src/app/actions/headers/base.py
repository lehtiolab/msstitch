from collections import OrderedDict


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
