import os
from app.readers import tsv as tsvreader
from app.dataformats import prottable as prottabledata


def create_lookup(fns, pqdb, poolnames, featcolnr, ms1_qcolpattern,
        isobqcolpattern, fdrcolpattern, flrcolpattern):
    psmnrpattern = prottabledata.HEADER_NO_PSMS_SUFFIX
    fullpsmpattern = prottabledata.HEADER_NO_FULLQ_PSMS
    patterns = [ms1_qcolpattern, fdrcolpattern, fullpsmpattern, flrcolpattern]
    storefuns = [pqdb.store_precursor_quants, pqdb.store_fdr, pqdb.store_fullq_psms,
            pqdb.store_ptm_flr]
    tablefn_map = create_tablefn_map(fns, pqdb, poolnames)
    feat_map = pqdb.get_feature_map()
    for pattern, storefun in zip(patterns, storefuns):
        if pattern is None:
            continue
        colmap = get_colmap(fns, pattern, single_col=True)
        if colmap:
            store_single_col_data(fns, tablefn_map, feat_map, storefun, colmap)
    if isobqcolpattern:
        isocolmap = get_colmap(fns, isobqcolpattern, antipattern=psmnrpattern)
        psmcolmap = get_colmap(fns, psmnrpattern)
        create_isobaric_quant_lookup(fns, tablefn_map, feat_map, pqdb, featcolnr,
                isocolmap, psmcolmap)


def create_tablefn_map(fns, pqdb, poolnames):
    """Stores protein/peptide table names in DB, returns a map with their
    respective DB IDs"""
    poolmap = {name: pid for (name, pid) in pqdb.get_all_poolnames()}
    pqdb.store_table_files([(poolmap[pool], os.path.basename(fn))
                            for fn, pool in zip(fns, poolnames)])
    return pqdb.get_tablefn_map()


def get_colmap(fns, pattern, single_col=False, antipattern=False):
    """For table files, loops through headers and checks which column(s)
    match a passed pattern. Those column(s) names are returned in a map with
    filenames as keys"""
    colmap = {}
    for fn in fns:
        header = tsvreader.get_tsv_header(fn)
        basefn = os.path.basename(fn)
        try:
            cols = tsvreader.get_cols_in_file(pattern, header, single_col)
        except RuntimeError:
            # Columns are not in this file
            cols = []
        if antipattern:
            try:
                anticols = tsvreader.get_cols_in_file(antipattern, header,
                                                      single_col)
            except RuntimeError:
                # The filtering "anti"-columns are not in the file, 
                anticols = []
            cols = [col for col in cols if col not in anticols]
        if cols:
            colmap[basefn] = cols
        else:
            return False
    return colmap


def store_single_col_data(fns, prottable_id_map, pacc_map, pqdbmethod, colmap):
    """General method to store single column data from protein tables
    in lookup"""
    to_store = []
    for fn, header, pquant in tsvreader.generate_tsv_pep_protein_quants(fns):
        pacc_id = pacc_map[pquant[header[0]]]
        pqdata = (pacc_id, prottable_id_map[fn], pquant[colmap[fn]])
        to_store.append(pqdata)
        if len(to_store) > 10000:
            pqdbmethod(to_store)
            to_store = []
    pqdbmethod(to_store)


def create_isobaric_quant_lookup(fns, tablefn_map, featmap, pqdb,
                                 featcolnr, allquantcols, psmcolmap):
    """Creates a lookup dict from peptide/protein quant input files and some
    input parameters. This assumes the order of quant columns and
    number-of-PSM columns is the same."""
    pqdb.store_quant_channels(
            map_psmnrcol_to_quantcol(allquantcols, psmcolmap, tablefn_map),
            psmcolmap)
    quantmap = pqdb.get_quantchannel_map()
    to_store = []
    for fn, header, pquant in tsvreader.generate_tsv_pep_protein_quants(fns):
        pqdata = get_isob_quant_data(pquant, header[featcolnr], tablefn_map[fn],
                                     featmap, quantmap)
        to_store.extend(pqdata)
        if len(to_store) > 10000:
            pqdb.store_isobaric_quants(to_store, psmcolmap)
            to_store = []
    pqdb.store_isobaric_quants(to_store, psmcolmap)


def map_psmnrcol_to_quantcol(quantcols, psmcols, tablefn_map):
    """This function yields tuples of table filename, isobaric quant column
    and if necessary number-of-PSM column"""
    if not psmcols:
        for fn in quantcols:
            for qcol in quantcols[fn]:
                yield (tablefn_map[fn], qcol)
    else:
        for fn in quantcols:
            for qcol, psmcol in zip(quantcols[fn], psmcols[fn]):
                yield (tablefn_map[fn], qcol, psmcol)


def get_isob_quant_data(pquant, featfield, fnid, featmap, qmap):
    """Turns a dict from a line of protein/peptide quant data into a set of
    tuples that can be stored in db"""
    feat_dbid = featmap[pquant[featfield]]
    for channel_name, (channel_id, psmfield) in qmap[fnid].items():
        if psmfield is None:
            yield (feat_dbid, channel_id, pquant[channel_name])
        else:
            yield (feat_dbid, channel_id, pquant[channel_name],
                   pquant[psmfield])


def build_proteintable(pqdb, mergecutoff, isobaric=False):
    """Fetches proteins and quants from joined lookup table, loops through
    them and when all of a protein's quants have been collected, yields the
    protein quant information."""
    # First get a dict with the data which is the same over sets (protein-gene mappings etc)
    pdmap = pqdb.create_pdata_map()
    previousfeat, outfeat = False, {}
    # Then loop over the SQLite output
    fdrfieldnr = pqdb.singlefields.index(prottabledata.HEADER_QVAL)
    # multi set output, means that a feature from DB can be delivered in
    # multiple lines (1/set). As soon as a new feat is received from db,
    # yield collected output and flush it for population with new protein
    for setname, feat_id, *setfeatvals in pqdb.merge_features():
        if feat_id != previousfeat and previousfeat and outfeat != {}:
            outfeat.update(pdmap[previousfeat])
            yield outfeat
            outfeat = {}
        if protein_pool_fdr_cutoff_ok(setfeatvals[fdrfieldnr], mergecutoff, feat_id, setname):
            for ix, field in enumerate(pqdb.singlefields):
                fieldname = '{}_{}'.format(setname, field)
                outfeat[fieldname] = setfeatvals[ix]
            if isobaric:
                channels = setfeatvals[len(pqdb.singlefields)].split(',')
                quants = setfeatvals[len(pqdb.singlefields)+1].split(',')
                try: 
                    ampsms = setfeatvals[len(pqdb.singlefields)+2].split(',')
                except AttributeError:
                    ampsms = [0] * len(channels)
                for ch, quant, ampsm in zip(channels, quants, ampsms):
                    outfeat['{}_{}'.format(setname, ch)] = quant
                    if ampsms:
                        outfeat['{}_{}{}'.format(setname, ch, prottabledata.HEADER_NO_PSMS_SUFFIX)] = ampsm
        previousfeat = feat_id
    if outfeat != {}:
        outfeat.update(pdmap[feat_id])
        yield outfeat



def protein_pool_fdr_cutoff_ok(fdr, fdrcutoff, feature, setname):
    if not fdrcutoff:
        return True
    warning = False
    if fdr is False:
        warning = True
    try:
        fdr = float(fdr)
    except (TypeError, ValueError):
        warning = True
    if warning:
        print('WARNING, filtering merge table on FDR but feat ID {} in '
                'set {} has FDR value "{}" in '
                'lookup.'.format(feature, setname, fdr))
        return False
    return fdr < fdrcutoff
