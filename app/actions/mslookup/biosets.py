def create_set_lookup(lookupdb, spectrafns, set_names):
    lookupdb.store_setnames(((x,) for x in set_names))
    set_id_map = lookupdb.get_setnames()
    biosets = ((fn, set_id_map[setname]) for fn, setname in
               zip(spectrafns, set_names))
    lookupdb.store_biosets(biosets)
    lookupdb.index_biosets()
