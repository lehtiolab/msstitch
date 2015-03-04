from decimal import Decimal, getcontext

from app.readers import spectra as specreader


def create_spectra_lookup(lookup, fn_spectra):
    """Stores all spectra rt and scan nr in db"""
    getcontext().prec = 14  # sets decimal point precision
    to_store = []
    mzmlmap = lookup.get_mzmlfile_map()
    for fn, spectrum, ns in fn_spectra:
        mzml_rt = float(Decimal(specreader.get_mzml_rt(spectrum, ns)))
        scan_nr = specreader.get_spec_scan_nr(spectrum)
        to_store.append((mzmlmap[fn], scan_nr, mzml_rt))
        if len(to_store) == 5000:
            lookup.store_mzmls(to_store)
            to_store = []
    lookup.store_mzmls(to_store)
    lookup.index_mzml()
