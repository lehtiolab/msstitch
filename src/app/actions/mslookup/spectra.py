DB_STORE_CHUNK = 500000


def create_spectra_lookup(lookup, fn_spectra):
    """Stores all spectra rt, injection time, and scan nr in db"""
    to_store = []
    mzmlmap = lookup.get_mzmlfile_map()
    for fn, spectrum in fn_spectra:
        spec_id = '{}_{}'.format(mzmlmap[fn], spectrum['scan'])
        mzml_rt = round(float(spectrum['rt']), 12)
        mzml_iit = round(float(spectrum['iit']), 12)
        mz = float(spectrum['mz'])
        to_store.append((spec_id, mzmlmap[fn], spectrum['scan'],
                         spectrum['charge'], mz, mzml_rt, mzml_iit))
        if len(to_store) == DB_STORE_CHUNK:
            lookup.store_mzmls(to_store)
            to_store = []
    lookup.store_mzmls(to_store)
    lookup.index_mzml()
