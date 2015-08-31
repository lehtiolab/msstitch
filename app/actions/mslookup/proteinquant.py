import re
import os
from app.readers import tsv as tsvreader


def create_proteinquant_lookup(fns, pqdb, poolnames, protacc_colnr,
                               ms1_qcolpattern=None, isobqcolpattern=None,
                               psmnrpattern=None, probcolpattern=None,
                               fdrcolpattern=None):
    poolmap = {name: pid for (name, pid) in pqdb.get_all_poolnames()}
    pqdb.store_protein_tables([(poolmap[pool], os.path.basename(fn))
                               for fn, pool in zip(fns, poolnames)])
    prottable_map = pqdb.get_protein_table_map()
    iso_quantcols, psmnrcolmap = {}, {}
    precur_quantcols, probcol, fdrcol = {}, {}, {}
    for fn in fns:
        header = tsvreader.get_tsv_header(fn)
        basefn = os.path.basename(fn)
        for colmap, pattern in zip([iso_quantcols, psmnrcolmap],
                                   [isobqcolpattern, psmnrpattern]):
            get_cols_in_file(colmap, pattern, basefn, header)
        for colmap, pattern in zip([precur_quantcols, probcol, fdrcol],
                                   [ms1_qcolpattern, probcolpattern,
                                    fdrcolpattern]):
            get_cols_in_file(colmap, pattern, basefn, header, single_col=True)
    if iso_quantcols and psmnrcolmap:
        create_isobaric_proteinquant_lookup(fns, prottable_map, pqdb,
                                            protacc_colnr,
                                            iso_quantcols, psmnrcolmap)
    if precur_quantcols:
        create_precursor_proteinquant_lookup(fns, prottable_map, pqdb,
                                             protacc_colnr, precur_quantcols)
    if probcol:
        create_probability_proteinquant_lookup(fns, prottable_map, pqdb,
                                               protacc_colnr, probcol)
    if fdrcol:
        create_fdr_proteinquant_lookup(fns, prottable_map, pqdb,
                                       protacc_colnr, fdrcol)


def create_protein_lookup(fns, prottable_id_map, pqdbmethod, protacc_colnr,
                          colmap):
    """General method to store single column data from protein tables
    in lookup"""
    to_store = []
    for fn, header, pquant in tsvreader.generate_tsv_protein_quants(fns):
        pqdata = (pquant[header[protacc_colnr]], prottable_id_map[fn],
                  pquant[colmap[fn]])
        to_store.append(pqdata)
        if len(to_store) > 10000:
            pqdbmethod(to_store)
            to_store = []
    pqdbmethod(to_store)


def create_probability_proteinquant_lookup(fns, prottable_map, pqdb,
                                           protacc_colnr, probcolmap):
    """Stores protein probability"""
    create_protein_lookup(fns, prottable_map, pqdb.store_protprob,
                          protacc_colnr, probcolmap)


def create_fdr_proteinquant_lookup(fns, prottable_map, pqdb,
                                   protacc_colnr, fdrcolmap):
    """Stores protein FDR"""
    create_protein_lookup(fns, prottable_map, pqdb.store_protfdr,
                          protacc_colnr, fdrcolmap)


def create_precursor_proteinquant_lookup(fns, prottable_map, pqdb,
                                         protacc_colnr, quantcolmap):
    """Stores protein precursor quant data"""
    create_protein_lookup(fns, prottable_map, pqdb.store_precursor_protquants,
                          protacc_colnr, quantcolmap)


def create_isobaric_proteinquant_lookup(fns, prottable_map, pqdb,
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
                                         protacc_colnr, quantmap)
        to_store.extend(pqdata)
        if len(to_store) > 10000:
            pqdb.store_isobaric_protquants(to_store)
            to_store = []
    pqdb.store_isobaric_protquants(to_store)


def get_cols_in_file(column_map, pattern, fn, header, single_col=False):
    if pattern is None:
        return {}
    cols_found = get_columns_by_pattern(header, pattern)
    if single_col:
        cols_found = cols_found[0]
    column_map.update({fn: cols_found})
    return column_map


def get_columns_by_pattern(header, pattern):
    columns = []
    for field in header:
        if re.search(pattern, field) is not None:
            columns.append(field)
    if not columns:
        raise RuntimeError('Could not find fieldname in header with '
                           'pattern: {}'.format(pattern))
    return columns


def map_psmnrcol_to_quantcol(quantcols, psmcols, prottable_map):
    for fn in quantcols:
        for qcol, psmcol in zip(quantcols[fn], psmcols[fn]):
            yield (prottable_map[fn], qcol, psmcol)


def get_isob_protquant_data(pquant, header, fnid, acccol, qmap):
    # (protein_acc, quantmap[qcol], quantvalue, amount_peptides)
    """Turns a dict from a line of protein quant data into a set of
    tuples that can be stored in db"""
    protacc = pquant[header[acccol]]
    for channel_name, (channel_id, psmfield) in qmap[fnid].items():
        yield (protacc, channel_id, pquant[channel_name], pquant[psmfield])
