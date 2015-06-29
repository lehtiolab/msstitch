from app.lookups.sqlite import quant as sqlite
from app.readers import tsv as readers
from app.dataformats import mzidtsv as mzidtsvdata


def generate_psms_quanted(quantdb, tsvfn, isob_header, oldheader,
                          is_ibariq=False, precursor=False): 
    """Takes dbfn and connects, gets quants for each line in tsvfn, sorts
    them in line by using keys in quantheader list."""
    allquants = quantdb.select_all_psm_quants()
    quant = next(allquants)
    for rownr, psm in enumerate(readers.generate_tsv_psms(tsvfn, oldheader)):
        outpsm = {x: y for x, y in psm.items()}
        if quant[3] is not None:
            outpsm.update({mzidtsvdata.HEADER_PRECURSOR_QUANT: str(quant[3])})
        else:
            outpsm.update({mzidtsvdata.HEADER_PRECURSOR_QUANT: 'NA'})
        isoquants = {}
        while quant[0] == rownr:
             isoquants.update({quant[1]: str(quant[2])})
             quant = next(allquants)
        outpsm.update(get_quant_NAs(isoquants, isob_header))
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


def create_tsv_header_quant(tsvfn, quantheader):
    """Returns tsvheader split list with quant header appended"""
    return readers.get_tsv_header(tsvfn) + [str(x[0]) for x in quantheader]


def get_quant_NAs(quantdata, quantheader):
    """Takes quantdata in a dict and header with quantkeys
    (eg iTRAQ isotopes). Returns dict of quant intensities
    with missing keys set to NA."""
    out = {}
    for qkey in quantheader:
        out[qkey] = quantdata.get(qkey, 'NA')
    return out
