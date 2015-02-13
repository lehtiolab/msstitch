from app.lookups.sqlite import quant as sqlite
from app.readers import tsv as readers
from app.dataformats import mzidtsv as mzidtsvdata


def generate_psms_quanted(quantdbfn, tsvfn, isob_header, oldheader,
                          is_ibariq=None, rttolerance=None, mztolerance=None):
    """Takes dbfn and connects, gets quants for each line in tsvfn, sorts
    them in line by using keys in quantheader list."""
    quantdb = sqlite.QuantDB(quantdbfn)
    quantfunctions = []
    if is_ibariq is not None:
        quantfunctions.append(lookup_iso_quant)
    if None not in [rttolerance, mztolerance]:
        quantfunctions.append(lookup_precursor_quant)
    for psm in readers.generate_tsv_psms(tsvfn, oldheader):
        outpsm = {x: y for x, y in psm.items()}
        specfile = outpsm[mzidtsvdata.HEADER_SPECFILE]
        scannr = outpsm[mzidtsvdata.HEADER_SCANNR]
        charge = outpsm[mzidtsvdata.HEADER_CHARGE]
        mz = outpsm[mzidtsvdata.HEADER_PRECURSOR_MZ]
        outpsm.update(lookup_quant(specfile, scannr, charge,
                      quantfunctions, mz, rttolerance, mztolerance,
                      quantdb, isob_header))
        yield outpsm


def get_full_and_isobaric_headers(oldheader, quantdbfn, isobaric=False,
                                  precursor=False):
    fullheader = oldheader
    if precursor:
        fullheader += mzidtsvdata.HEADER_PRECURSOR_QUANT
    if isobaric:
        quantdb = sqlite.QuantDB(quantdbfn)
        quantmap = quantdb.get_all_quantmaps()
        isob_header = sorted([x[0] for x in quantmap])
        fullheader += isob_header
    else:
        isob_header = None
    return fullheader, isob_header


def get_quant_header(oldheader, quantdbfn):
    quantdb = sqlite.QuantDB(quantdbfn)
    quantmap = quantdb.get_all_quantmaps()
    return oldheader + sorted([x[0] for x in quantmap])


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
                 rttol, mztol, quantdb, isob_header=None):
    outquants = {}
    for func in quantfunctions:
        outquants.update(func(quantdb, specfile, scannr, charge,
                              mz, rttol, mztol, header=isob_header))


def lookup_iso_quant(quantdb, spectrafile, scannr, *args, **kwargs):
    """Outputs dict with keys == quantname, values == quantintensity."""
    dbquants = quantdb.lookup_isobaric_quant(spectrafile, scannr)
    return get_quant_NAs({x[0]: str(x[1]) for x in dbquants}, kwargs['header'])


def lookup_precursor_quant(quantheader, quantdb, spectrafile, scannr,
                           charge, mz, rttol, mztol, **kwargs):
    """Lookup quant features in db that lie inside the m/z and retention time
    tolerance limits. Returns the one feature which has the best matching
    m/z, but not necessarily best matching retention time"""
    def get_minmax(center, tolerance):
        center, tolerance = float(center), float(tolerance)
        return center - tolerance, center + tolerance
    ms2_rt = quantdb.lookup_retention_time(spectrafile, scannr)
    minrt, maxrt = get_minmax(ms2_rt, rttol)
    minmz, maxmz = get_minmax(mz, mztol)
    dbquants = quantdb.lookup_precursor_quant(spectrafile, charge, minrt,
                                              maxrt, minmz, maxmz)
    # m/z has index 1 from db output tuple
    features = {abs(mz - x[1]): x for x in dbquants}
    return {mzidtsvdata.HEADER_PRECURSOR_QUANT: features[max(features)]}
