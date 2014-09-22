from app import sqlite
from app.readers import tsv as readers


def generate_psms_quanted(quantdbfn, tsvfn, quantheader, oldheader):
    """Takes dbfn and connects, gets quants for each line in tsvfn, sorts
    them in line by using keys in quantheader list."""
    quantdb = sqlite.QuantDB(quantdbfn)
    for line, (specfile,
               scannr) in readers.get_mzidtsv_lines_scannr_specfn(tsvfn):
        outlinedict = {k: v for k, v in zip(oldheader, line)}
        quantdata = lookup_quant(specfile, scannr, quantdb)
        quantlinedict = get_quant_NAs(quantdata, quantheader)
        outlinedict.update(quantlinedict)
        yield outlinedict


def get_quant_header(oldheader, quantdbfn):
    quantdb = sqlite.QuantDB(quantdbfn)
    quantmap = quantdb.get_all_quantmaps()
    return oldheader + sorted([x[0] for x in quantmap])


def create_tsv_header_quant(tsvfn, quantheader):
    """Returns tsvheader split list with quant header appended"""
    return readers.get_tsv_header(tsvfn) + [str(x[0]) for x in quantheader]


def get_quant_NAs(quantdata, quantheader):
    """Takes quantdata in a dict and header with quantkeys
    (eg iTRAQ isotopes). Returns dict of quant intensities
    with missing keys set to NA."""
    out = {}
    for qkey in quantheader:
        out[qkey] = quantdata.get(qkey, 'NA')
    return out


def lookup_quant(spectrafile, psm_scan_nr, quantdb):
    """Outputs dict with keys == quantname, values == quantintensity."""
    dbquants = quantdb.lookup_quant(spectrafile, psm_scan_nr)
    return {x[0]: x[1] for x in dbquants}
