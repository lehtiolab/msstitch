from app import sqlite
from app import formatting
from app.readers import ompstitch as readers


def generate_psms_quanted(quantdbfn, psms):
    quantdb = sqlite.QuantDB(quantdbfn)
    for psm in psms:
        psm_scan_nr = readers.get_psm_scan_nr(psm)
        quant = lookup_quant(psm_scan_nr, quantdb)
        output = create_psm_out(psm, quant)
        formatting.clear_el(psm)
        yield output


def lookup_quant(psm_scan_nr, quantdb):
    dbquants = quantdb.find_quants(psm_scan_nr)
    return {x[0]: x[1] for x in dbquants}


def create_psm_out(psm, quant):
    """Return dict with keys == headerfields,
    values == psm values for those."""
    # FIXME look how Zazzi's mzid with Percolator turns out
    psmout = {}
    header = ['mzIdentML ID', 'scan_nr', 'mzml_filename', 'sequence',
              'q-value', 'PEP']  # protein IDs, peptide ID in percolator out?
    return header, psmout
