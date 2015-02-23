from app.lookups.sqlite import quant as sqlite
from app.readers import tsv as readers
from app.dataformats import mzidtsv as mzidtsvdata


def generate_psms_quanted(quantdbfn, tsvfn, isob_header, oldheader,
                          is_ibariq=False, rttolerance=None, mztolerance=None,
                          mztoltype=None, spec_column=None):
    """Takes dbfn and connects, gets quants for each line in tsvfn, sorts
    them in line by using keys in quantheader list."""
    quantdb = sqlite.QuantDB(quantdbfn)
    quantfunctions = []
    if is_ibariq:
        quantfunctions.append(lookup_iso_quant)
    if None not in [rttolerance, mztolerance, mztoltype]:
        quantfunctions.append(lookup_precursor_quant)
    for psm in readers.generate_tsv_psms(tsvfn, oldheader):
        outpsm = {x: y for x, y in psm.items()}
        if spec_column is not None:
            specfile = outpsm[oldheader[spec_column - 1]]
        else:
            specfile = outpsm[mzidtsvdata.HEADER_SPECFILE]
        scannr = outpsm[mzidtsvdata.HEADER_SCANNR]
        charge = outpsm[mzidtsvdata.HEADER_CHARGE]
        mz = outpsm[mzidtsvdata.HEADER_PRECURSOR_MZ]
        outpsm.update(lookup_quant(specfile, scannr, charge,
                      quantfunctions, mz, rttolerance, mztolerance,
                      mztoltype, quantdb, isob_header))
        yield outpsm


def get_full_and_isobaric_headers(oldheader, quantdbfn, isobaric=False,
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
    # driven stuff. then output here oldheader and new fields as tuple.
    fullheader = oldheader
    if precursor:
        fullheader += [mzidtsvdata.HEADER_PRECURSOR_QUANT]
    if isobaric:
        quantdb = sqlite.QuantDB(quantdbfn)
        quantmap = quantdb.get_all_quantmaps()
        isob_header = sorted([x[0] for x in quantmap])
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


def lookup_quant(specfile, scannr, charge, quantfunctions, mz,
                 rttol, mztol, mztoltype, quantdb, isob_header=None):
    outquants = {}
    for func in quantfunctions:
        outquants.update(func(quantdb, specfile, scannr, charge,
                              mz, rttol, mztol, mztoltype, header=isob_header))
    return outquants


def lookup_iso_quant(quantdb, spectrafile, scannr, *args, **kwargs):
    """Outputs dict with keys == quantname, values == quantintensity."""
    dbquants = quantdb.lookup_isobaric_quant(spectrafile, scannr)
    return get_quant_NAs({x[0]: str(x[1]) for x in dbquants}, kwargs['header'])


def lookup_precursor_quant(quantdb, spectrafile, scannr,
                           charge, mz, rttol, mztol, mztoltype, **kwargs):
    """Lookup quant features in db that lie inside the m/z and retention time
    tolerance limits. Returns the one feature which has the best matching
    m/z, but not necessarily best matching retention time"""
    def get_minmax(center, tolerance, toltype=None):
        center = float(center)
        if toltype == 'ppm':
            tolerance = int(tolerance) / 1000000 * center
        elif toltype == 'Da':
            tolerance = float(tolerance)
        return center - tolerance, center + tolerance
    ms2_rt = quantdb.lookup_retention_time(spectrafile, scannr)[0][0]
    minrt, maxrt = get_minmax(ms2_rt, rttol / 60)
    minmz, maxmz = get_minmax(mz, mztol, mztoltype)
    dbquants = quantdb.lookup_precursor_quant(spectrafile, charge, minrt,
                                              maxrt, minmz, maxmz)
    # m/z has index 0 from db output tuple
    features = {abs(float(mz) - x[0]): x[1] for x in dbquants}
    if features:
        outquant = str(features[min(features)])
    else:
        outquant = 'NA'
    return {mzidtsvdata.HEADER_PRECURSOR_QUANT: outquant}
