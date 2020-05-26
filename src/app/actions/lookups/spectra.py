import os

DB_STORE_CHUNK = 500000


def create_bioset_lookup(lookupdb, spectrafns, set_names):
    """Fills lookup database with biological set names"""
    unique_setnames = set(set_names)
    lookupdb.store_biosets(((x,) for x in unique_setnames))
    set_id_map = lookupdb.get_setnames()
    mzmlfiles = ((os.path.basename(fn), set_id_map[setname])
                 for fn, setname in zip(spectrafns, set_names))
    lookupdb.store_mzmlfiles(mzmlfiles)
    lookupdb.index_biosets()


def create_spectra_lookup(lookup, fn_spectra):
    """Stores all spectra rt, injection time, and scan ID in db"""
    dbspectra, dbioninj, dbionmob = [], [], []
    count = 0
    mzmlmap = lookup.get_mzmlfile_map()
    for fn, spectrum in fn_spectra:
        spec_id = '{}_{}'.format(mzmlmap[fn], spectrum['specscanid'])
        mzml_rt = round(float(spectrum['rt']), 12)
        mzml_iit = round(float(spectrum['iit']), 12)
        mzml_mob = round(float(spectrum['ionmob']), 12)
        mz = float(spectrum['mz'])
        dbspectra.append((spec_id, mzmlmap[fn], spectrum['specscanid'],
                         spectrum['charge'], mz, mzml_rt))
        dbioninj.append((spec_id, mzml_iit)) if mzml_iit else False
        dbionmob.append((spec_id, mzml_mob)) if mzml_mob else False
        count += 1
        if len(dbspectra) == DB_STORE_CHUNK:
            lookup.store_mzmls(dbspectra, dbioninj, dbionmob)
            dbspectra, dbioninj, dbionmob  = [], [], []
    lookup.store_mzmls(dbspectra, dbioninj, dbionmob)
    lookup.index_mzml()
    print('Stored {} spectra in lookup'.format(count))
