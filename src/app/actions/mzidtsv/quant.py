from app.readers import tsv as readers
from app.dataformats import mzidtsv as mzidtsvdata


def generate_psms_quanted(quantdb, tsvfn, isob_header, oldheader,
                          isobaric=False, precursor=False):
    """Takes dbfn and connects, gets quants for each line in tsvfn, sorts
    them in line by using keys in quantheader list."""
    allquants, sqlfields = quantdb.select_all_psm_quants(isobaric, precursor)
    quant = next(allquants)
    for rownr, psm in enumerate(readers.generate_tsv_psms(tsvfn, oldheader)):
        outpsm = {x: y for x, y in psm.items()}
        if precursor:
            pquant = quant[sqlfields['precursor']]
            if pquant is None:
                pquant = 'NA'
            outpsm.update({mzidtsvdata.HEADER_PRECURSOR_QUANT: str(pquant)})
        if isobaric:
            isoquants = {}
            while quant[0] == rownr:
                isoquants.update({quant[sqlfields['isochan']]:
                                  str(quant[sqlfields['isoquant']])})
                try:
                    quant = next(allquants)
                except StopIteration:
                    # last PSM, break from while loop or it is not yielded at all
                    break
            outpsm.update(get_quant_NAs(isoquants, isob_header))
        else:
            try:
                quant = next(allquants)
            except StopIteration:
                # last PSM, needs explicit yield/break or it will not be yielded
                yield outpsm
                break
        yield outpsm


def get_full_and_isobaric_headers(oldheader, quantdb, isobaric=False,
                                  precursor=False):
    # FIXME:
    # if we let db decide, and already isobaric done on tsv gets a double db,
    # then automatically we get twice isobaric on the tsv header. not good,
    # because it will fuck up the psm dicts per line.
    # also, maybe we should make sure there is no overwriting the header fields
    # with new fields in case that would  happne, or output a set.
    # is there any other scenario where we dont want a specific part of quant
    # data included in the tsv except 'it is already there'?
    # is we're outputting a set, we should do this as a general method for tsv
    # driven stuff. then output here oldheader and new fields as tuple
    # FIXME Make sure header is in mass-order, not alphabetical as it is now
    fullheader = oldheader
    if precursor:
        fullheader += [mzidtsvdata.HEADER_PRECURSOR_QUANT]
    if isobaric:
        isob_header = [x[0] for x in quantdb.get_all_quantmaps()]
        fullheader += isob_header
    else:
        isob_header = None
    return fullheader, isob_header


def get_quant_NAs(quantdata, quantheader):
    """Takes quantdata in a dict and header with quantkeys
    (eg iTRAQ isotopes). Returns dict of quant intensities
    with missing keys set to NA."""
    out = {}
    for qkey in quantheader:
        out[qkey] = quantdata.get(qkey, 'NA')
    return out
