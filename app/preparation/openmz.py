from app import sqlite
from app.readers import tsv as readers


def generate_psms_quanted(quantdbfn, tsvfn, quantheader):
    """Takes dbfn and connects, gets quants for each line in tsvfn, sorts
    them in line by using keys in quantheader list."""
    quantdb = sqlite.QuantDB(quantdbfn)
    for line, (specfile,
               scannr) in readers.get_mzidtsv_lines_scannr_specfn(tsvfn):
        quantdata = lookup_quant(specfile, scannr, quantdb)
        quantline = convert_quantdata_to_line(quantdata, quantheader)
        yield line + quantline


def get_quant_header(quantdbfn):
    quantdb = sqlite.QuantDB(quantdbfn)
    quantmap = quantdb.get_all_quantmaps()
    return sorted(quantmap)


def create_tsv_header_quant(tsvfn, quantheader):
    """Returns tsvheader split list with quant header appended"""
    return next(readers.get_tsv_header(tsvfn)) + quantheader


def convert_quantdata_to_line(quantdata, quantheader):
    """Takes quantdata in a dict and header with quantkeys
    (eg iTRAQ isotopes). Returns list of quant intensities from dict.
    NA if quantheader item was not a key in the quantdata dict."""
    line = []
    for qkey in quantheader:
        try:
            line.append(quantdata[qkey])
        except KeyError:
            line.append('NA')
    return line


def lookup_quant(spectrafile, psm_scan_nr, quantdb):
    """Outputs dict with keys == quantname, values == quantintensity."""
    dbquants = quantdb.lookup_quant(spectrafile, psm_scan_nr)
    return {x[0]: x[1] for x in dbquants}



