import re
import os
from app.readers import tsv as tsvreader


def create_proteinquant_lookup(fns, pqdb, protacc_colnr,
                               ms1_qcolpattern=None, isobqcolpattern=None,
                               psmnrpattern=None):
    iso_quantcols = {}
    precur_quantcol = {}
    psmnrcolmap = {}
    for fn in fns:
        header = tsvreader.get_tsv_header(fn)
        basefn = os.path.basename(fn)
        iso_quantcols.update({basefn: get_quantcolumns(header,
                                                       isobqcolpattern)})
        psmnrcolmap.update({basefn: get_quantcolumns(header, psmnrpattern)})
        precur_quantcol.update({basefn: get_quantcolumns(header,
                                                         ms1_qcolpattern)[0]})
    if isobqcolpattern is not None and psmnrpattern is not None:
        create_isobaric_proteinquant_lookup(fns, pqdb, protacc_colnr,
                                            iso_quantcols, psmnrcolmap)
    if ms1_qcolpattern is not None:
        create_precursor_proteinquant_lookup(fns, pqdb, protacc_colnr,
                                             precur_quantcol)


def create_precursor_proteinquant_lookup(fns, pqdb, protacc_colnr,
                                         quantcolmap):
    to_store = []
    for fn, header, pquant in tsvreader.generate_tsv_protein_quants(fns):
        pqdata = (pquant[header[protacc_colnr]], pquant[quantcolmap[fn]])
        to_store.append(pqdata)
        if len(to_store) > 10000:
            pqdb.store_precursor_protquants(to_store)
            to_store = []
    pqdb.store_precursor_protquants(to_store)


def create_isobaric_proteinquant_lookup(fns, pqdb, protacc_colnr,
                                        allquantcols, psmnrcolmap):
    """Creates a lookup dict from protein quant input files and some
    input parameters. This assumes the order of quant columns and
    number-of-PSM columns is the same."""
    pqdb.store_quant_channels(map_psmnrcol_to_quantcol(allquantcols,
                                                       psmnrcolmap))
    quantmap = pqdb.get_quantchannel_map()
    to_store = []
    for fn, header, pquant in tsvreader.generate_tsv_protein_quants(fns):
        pqdata = get_isob_protquant_data(pquant, header, fn, protacc_colnr,
                                         quantmap)
        to_store.extend(pqdata)
        if len(to_store) > 10000:
            pqdb.store_isobaric_protquants(to_store)
            to_store = []
    pqdb.store_isobaric_protquants(to_store)


def get_quantcolumns(header, quantcolpattern):
    if quantcolpattern is None:
        return None
    columns = []
    for field in header:
        if re.search(quantcolpattern, field) is not None:
            columns.append(field)
    return columns


def map_psmnrcol_to_quantcol(quantcols, psmcols):
    for fn in quantcols:
        for qcol, psmcol in zip(quantcols[fn], psmcols[fn]):
            yield (fn, qcol, psmcol)


def get_isob_protquant_data(pquant, header, fn, acccol, qmap):
    # (protein_acc, quantmap[qcol], quantvalue, amount_peptides)
    """Turns a dict from a line of protein quant data into a set of
    tuples that can be stored in db"""
    protacc = pquant[header[acccol]]
    for channel_name, (channel_id, psmfield) in qmap[fn].items():
        yield (protacc, channel_id, pquant[channel_name], pquant[psmfield])
