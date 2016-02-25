import os
from app.readers import tsv as tsvreader


def create_peptidequant_lookup(fns, pqdb, poolnames, pepseq_colnr,
                               ms1_qcolpattern=None, isobqcolpattern=None,
                               psmnrpattern=None, fdrcolpattern=None,
                               pepcolpattern=None):
    """Calls lower level function to create a peptide quant lookup"""
    patterns = [ms1_qcolpattern, fdrcolpattern, pepcolpattern]
    storefuns = [pqdb.store_precursor_quants, pqdb.store_fdr,
                 pqdb.store_pep]
    create_pep_protein_quant_lookup(fns, pqdb, poolnames, pepseq_colnr,
                                    patterns, storefuns,
                                    isobqcolpattern, psmnrpattern)


def create_proteinquant_lookup(fns, pqdb, poolnames, protacc_colnr,
                               ms1_qcolpattern=None, isobqcolpattern=None,
                               psmnrpattern=None, probcolpattern=None,
                               fdrcolpattern=None, pepcolpattern=None):
    """Calls lower level function to create a protein quant lookup"""
    patterns = [ms1_qcolpattern, probcolpattern, fdrcolpattern, pepcolpattern]
    storefuns = [pqdb.store_precursor_quants, pqdb.store_probability,
                 pqdb.store_fdr, pqdb.store_pep]
    create_pep_protein_quant_lookup(fns, pqdb, poolnames, protacc_colnr,
                                    patterns, storefuns, isobqcolpattern,
                                    psmnrpattern)


def create_pep_protein_quant_lookup(fns, pqdb, poolnames, featcolnr, patterns,
                                    storefuns, isobqcolpattern=None,
                                    psmnrpattern=None):
    """Does the work when creating peptide and protein quant lookups. This
    loops through storing options and parses columns, passing on to the
    storing functions"""
    tablefn_map = create_tablefn_map(fns, pqdb, poolnames)
    feat_map = pqdb.get_feature_map()
    for pattern, storefun in zip(patterns, storefuns):
        if pattern is None:
            continue
        colmap = get_colmap(fns, pattern, single_col=True)
        if colmap:
            store_single_col_data(fns, tablefn_map, feat_map,
                                  storefun, featcolnr, colmap)
    if isobqcolpattern is not None:
        isocolmap = get_colmap(fns, isobqcolpattern, antipattern=psmnrpattern)
    else:
        return
    if psmnrpattern is not None:
        psmcolmap = get_colmap(fns, psmnrpattern)
    create_isobaric_quant_lookup(fns, tablefn_map,
                                 feat_map, pqdb,
                                 featcolnr,
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
        cols = tsvreader.get_cols_in_file(pattern, header, single_col)
        if antipattern:
            anticols = tsvreader.get_cols_in_file(antipattern, header,
                                                  single_col)
            cols = [col for col in cols if col not in anticols]
        if cols:
            colmap[basefn] = cols
    return colmap


def store_single_col_data(fns, prottable_id_map, pacc_map, pqdbmethod,
                          protacc_colnr, colmap):
    """General method to store single column data from protein tables
    in lookup"""
    to_store = []
    for fn, header, pquant in tsvreader.generate_tsv_pep_protein_quants(fns):
        pacc_id = pacc_map[pquant[header[protacc_colnr]]]
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
    pqdb.store_quant_channels(map_psmnrcol_to_quantcol(allquantcols,
                                                       psmcolmap,
                                                       tablefn_map))
    quantmap = pqdb.get_quantchannel_map()
    to_store = []
    for fn, header, pquant in tsvreader.generate_tsv_pep_protein_quants(fns):
        pqdata = get_isob_quant_data(pquant, header, tablefn_map[fn],
                                     featmap, featcolnr, quantmap)
        to_store.extend(pqdata)
        if len(to_store) > 10000:
            pqdb.store_isobaric_quants(to_store)
            to_store = []
    pqdb.store_isobaric_quants(to_store)


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


def get_isob_quant_data(pquant, header, fnid, featmap, featcol, qmap):
    """Turns a dict from a line of protein/peptide quant data into a set of
    tuples that can be stored in db"""
    feat_dbid = featmap[pquant[header[featcol]]]
    for channel_name, (channel_id, psmfield) in qmap[fnid].items():
        if psmfield is None:
            yield (feat_dbid, channel_id, pquant[channel_name])
        else:
            yield (feat_dbid, channel_id, pquant[channel_name],
                   pquant[psmfield])
