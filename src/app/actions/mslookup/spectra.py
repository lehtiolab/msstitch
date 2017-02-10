from decimal import Decimal, getcontext

DB_STORE_CHUNK = 500000


def create_spectra_lookup(lookup, fn_spectra):
    """Stores all spectra rt and scan nr in db"""
    getcontext().prec = 14  # sets decimal point precision
    to_store = []
    mzmlmap = lookup.get_mzmlfile_map()
    for fn, spectrum in fn_spectra:
        spec_id = '{}_{}'.format(mzmlmap[fn], spectrum['scan'])
        mzml_rt = float(Decimal(spectrum['rt']))
        mz = float(Decimal(spectrum['mz']))
        to_store.append((spec_id, mzmlmap[fn], spectrum['scan'],
                         spectrum['charge'], mz, mzml_rt))
        if len(to_store) == DB_STORE_CHUNK:
            lookup.store_mzmls(to_store)
            to_store = []
    lookup.store_mzmls(to_store)
    lookup.index_mzml()
