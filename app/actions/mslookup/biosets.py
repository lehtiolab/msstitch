import os


def create_bioset_lookup(lookupdb, spectrafns, set_names):
    lookupdb.store_biosets(((x,) for x in set_names))
    set_id_map = lookupdb.get_setnames()
    mzmlfiles = ((os.path.basename(fn), set_id_map[setname])
                 for fn, setname in zip(spectrafns, set_names))
    lookupdb.store_mzmlfiles(mzmlfiles)
    lookupdb.index_biosets()
