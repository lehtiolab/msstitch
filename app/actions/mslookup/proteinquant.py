import re
from app.readers import tsv as tsvreader


def create_proteinquant_lookup(fns, pqdb, protacc_colnr, qcolpattern, psmnrpattern):
    """Creates a lookup dict from protein quant input files and some
    input parameters. This assumes the order of quant columns and 
    number-of-PSM columns is the same."""
    allquantcols = {}
    psmnrcolmap = {}
    for fn in fns:
        header = tsvreader.get_tsv_header(fn)
        allquantcols.update({fn: get_quantcolumns(header, qcolpattern)})
        psmnrcolmap.update({fn: get_quantcolumns(header, psmnrpattern)})
    pqdb.store_quant_channels(map_psmnrcol_to_quantcol(allquantcols,
                                                       psmnrcolmap))
    quantmap = pqdb.get_quantchannel_map()
    to_store = []
    for fn, header, pquant in tsvreader.generate_tsv_protein_quants(fns):
        pqdata = get_protquant_data(pquant, header, fn, protacc_colnr, 
                                    quantmap)
        to_store.extend(pqdata)
        if len(to_store) > 10000:
            pqdb.store_protquants(to_store)
            to_store = []
    pqdb.store_protquants(to_store)


def get_quantcolumns(header, quantcolpattern):
    columns = []
    for field in header:
        if re.search(quantcolpattern, field) is not None:
            columns.append(field)
    return columns


def map_psmnrcol_to_quantcol(quantcols, psmcols):
    for fn in quantcols:
        for qcol, psmcol in zip(quantcols[fn], psmcols[fn]):
            yield (fn, qcol, psmcol)


def get_protquant_data(pquant, header, fn, acccol, qmap):
    # (protein_acc, quantmap[qcol], quantvalue, amount_peptides)
    """Turns a dict from a line of protein quant data into a set of
    tuples that can be stored in db"""
    protacc = pquant[header[acccol]]
    for channel_name, (channel_id, psmfield) in qmap[fn].items():
        yield (protacc, channel_id, pquant[channel_name], pquant[psmfield])
