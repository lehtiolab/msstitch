from app.readers import tsv as tsvreader


def create_proteinquant_lookup(fns, pqdb, protacc_colnr, qcolpattern):
    """Creates a lookup dict from protein quant input files and some
    input parameters"""
    allquantcols = set()
    for fn in fns:
        header = tsvreader.get_tsv_header(fn)
        allquantcols.update([header[x] for x in
                             get_quantcolumns(header, qcolpattern)])
    pqdb.store_quant_channels(allquantcols)
    quantmap = pqdb.get_quantchannel_ids()
    to_store = []
    for pquant in tsvreader.generate_protein_quants(fns):
        pqdata = get_protquant_data(pquant, header, protacc_colnr, qcolpattern,
                                    quantmap)
        to_store.extend(pqdata)
        if len(to_store) > 5000:
            pqdb.store_protquants(to_store)
            to_store = []
    pqdb.store_protquants(to_store)


def get_quantcolumns(header, quantcolpattern):
    for field in header:
        if quantcolpattern in field:
            yield field


def get_protquant_data(pquant, header, acccol, qcolpattern, qmap):
    # (protein_acc, quantmap[qcol], quantvalue, amount_peptides)
    """Turns a dict from a line of protein quant data into a set of
    tuples that can be stored in db"""
    protacc = pquant[header[acccol]]
    qfields = get_quantcolumns(header, qcolpattern)
    for field in qfields:
        yield (protacc, qmap[field], pquant[field])
