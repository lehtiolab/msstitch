import os
from app.readers import tsv as tsvreader


def create_tablefn_map(fns, pqdb, poolnames):
    poolmap = {name: pid for (name, pid) in pqdb.get_all_poolnames()}
    pqdb.store_table_files([(poolmap[pool], os.path.basename(fn))
                            for fn, pool in zip(fns, poolnames)])
    return pqdb.get_tablefn_map()


def get_colmap(fns, pattern, single_col=False):
    colmap = {}
    for fn in fns:
        header = tsvreader.get_tsv_header(fn)
        basefn = os.path.basename(fn)
        cols = tsvreader.get_cols_in_file(pattern, header, single_col)
        if cols:
            colmap[basefn] = cols
    return colmap


def create_proteinquant_lookup(fns, pqdb, poolnames, protacc_colnr,
                               ms1_qcolpattern=None, isobqcolpattern=None,
                               psmnrpattern=None, probcolpattern=None,
                               fdrcolpattern=None, pepcolpattern=None):
    prottable_map = create_tablefn_map(fns, pqdb, poolnames)
    protein_acc_map = pqdb.get_protein_acc_map()
    patterns = [ms1_qcolpattern, probcolpattern, fdrcolpattern, pepcolpattern]
    storefuns = [pqdb.store_precursor_protquants, pqdb.store_protprob,
                 pqdb.store_protfdr, pqdb.store_protpep]
    for pattern, storefun in zip(patterns, storefuns):
        if pattern is None:
            continue
        colmap = get_colmap(fns, pattern, single_col=True)
        if colmap:
            create_protein_lookup(fns, prottable_map, protein_acc_map,
                                  storefun, protacc_colnr, colmap)
    if isobqcolpattern is not None and psmnrpattern is not None:
        isocolmap = get_colmap(fns, isobqcolpattern)
        psmcolmap = get_colmap(fns, psmnrpattern)
        if isocolmap and psmcolmap:
            create_isobaric_proteinquant_lookup(fns, prottable_map,
                                                protein_acc_map, pqdb,
                                                protacc_colnr,
                                                isocolmap, psmcolmap)


def create_protein_lookup(fns, prottable_id_map, pacc_map, pqdbmethod,
                          protacc_colnr, colmap):
    """General method to store single column data from protein tables
    in lookup"""
    to_store = []
    for fn, header, pquant in tsvreader.generate_tsv_protein_quants(fns):
        pacc_id = pacc_map[pquant[header[protacc_colnr]]]
        pqdata = (pacc_id, prottable_id_map[fn], pquant[colmap[fn]])
        to_store.append(pqdata)
        if len(to_store) > 10000:
            pqdbmethod(to_store)
            to_store = []
    pqdbmethod(to_store)


def create_isobaric_proteinquant_lookup(fns, prottable_map, pacc_map, pqdb,
                                        protacc_colnr, allquantcols,
                                        psmnrcolmap):
    """Creates a lookup dict from protein quant input files and some
    input parameters. This assumes the order of quant columns and
    number-of-PSM columns is the same."""
    pqdb.store_quant_channels(map_psmnrcol_to_quantcol(allquantcols,
                                                       psmnrcolmap,
                                                       prottable_map))
    quantmap = pqdb.get_quantchannel_map()
    to_store = []
    for fn, header, pquant in tsvreader.generate_tsv_protein_quants(fns):
        pqdata = get_isob_protquant_data(pquant, header, prottable_map[fn],
                                         pacc_map, protacc_colnr, quantmap)
        to_store.extend(pqdata)
        if len(to_store) > 10000:
            pqdb.store_isobaric_protquants(to_store)
            to_store = []
    pqdb.store_isobaric_protquants(to_store)


def map_psmnrcol_to_quantcol(quantcols, psmcols, prottable_map):
    for fn in quantcols:
        for qcol, psmcol in zip(quantcols[fn], psmcols[fn]):
            yield (prottable_map[fn], qcol, psmcol)


def get_isob_protquant_data(pquant, header, fnid, pacc_map, acccol, qmap):
    # (protein_acc, quantmap[qcol], quantvalue, amount_peptides)
    """Turns a dict from a line of protein quant data into a set of
    tuples that can be stored in db"""
    protacc_id = pacc_map[pquant[header[acccol]]]
    for channel_name, (channel_id, psmfield) in qmap[fnid].items():
        yield (protacc_id, channel_id, pquant[channel_name], pquant[psmfield])
