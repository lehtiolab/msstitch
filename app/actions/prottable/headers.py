def get_header(oldheader=False, quant_psm_channels=False, addprotein_data=False,
               prottable_filenames=False, precursorarea=False, probability=False):
    if not oldheader:
        header = [prottabledata.HEADER_PROTEIN]
    else:
        header = oldheader[:]
    if quant_psm_channels:
        header.extend(quant_psm_channels)
    if addprotein_data:
        header = get_header_with_proteindata(header, prottable_filenames)
    if prottable_filenames or precursorarea:
        header = get_header_with_precursorarea(header, prottable_filenames)
    if probability:
        header = get_header_with_prot_probability(header, prottable_filenames)
    return header


def build_pool_header_field(pool, field):
    return '{}_{}'.format(pool, field)


def get_header_with_proteindata(header):
    ix = header.index(prottabledata.HEADER_PROTEIN) + 1
    new_data = [prottabledata.HEADER_DESCRIPTION,
                prottabledata.HEADER_COVERAGE,
                prottabledata.HEADER_NO_PROTEIN,
                prottabledata.HEADER_NO_UNIPEP,
                prottabledata.HEADER_NO_PEPTIDE,
                prottabledata.HEADER_NO_PSM,
                #prottabledata.HEADER_NO_QUANT_PSM,
                #prottabledata.HEADER_CV_QUANT_PSM,
                ]
    return header[:ix] + new_data + header[ix:]


def get_header_with_prot_probability(header, fns=False):
    ix = header.index(prottabledata.HEADER_PROTEIN) + 1
    if fns:
        new_data = [build_pool_header_field(fn, prottabledata.HEADER_PROBABILITY) for fn in fns]
    else:
        new_data = [prottabledata.HEADER_PROBABILITY]
    return header[:ix] + new_data + header[ix:]


def get_header_with_precursorarea(header, fns=False):
    ix = header.index(prottabledata.HEADER_PROTEIN) + 1
    if fns:
        quant_fields = [build_pool_header_field(fn, prottabledata.HEADER_AREA) for fn in fns]
    else:
        quant_fields = [prottabledata.HEADER_AREA]
    return header[:ix] + quant_fields + header[ix:]
