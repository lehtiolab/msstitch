from app import sqlite
from app import formatting
from app.readers import ompstitch as reader


def generate_psms_quanted(quantdbfn, psms):
    quantdb = sqlite.QuantDB(quantdbfn)
    for psm in psms:
        psm_scan_nr = reader.get_psm_scan_nr(psm)
        quant = lookup_quant(psm_scan_nr, quantdb)
        output = create_psm_out(psm, quant)
        formatting.clear_el(psm)
        yield output


def lookup_quant(psm_scan_nr, quantdb):
    dbquants = quantdb.find_quants(psm_scan_nr)
    return {x[0]: x[1] for x in dbquants}


def create_psm_out(psm, quant):
    # FIXME look how Zazzi's mzid with Percolator turns out
    pass
