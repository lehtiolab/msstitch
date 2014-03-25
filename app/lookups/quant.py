from app import sqlite
from app import formatting
from app.readers import ompstitch as readers


def create_quant_lookup(spectra, consensus_els):
    """Creates an sqlite lookup table of scannrs with quant data.

    spectra - an iterable of spectra
    consensus_els - a iterable with consensusElements"""
    quantdb = sqlite.QuantDB()
    quantdb.create_quantdb()
    quants = {}
    count = 0
    for consensus_el in consensus_els:
        count += 1
        rt = readers.get_consxml_rt(consensus_el)
        for spectrum in spectra:
            if rt == readers.get_mzml_rt(spectrum):
                spec_scan_nr = readers.get_spec_scan_nr(spectrum)
                break

        quants[spec_scan_nr] = get_quant_data(consensus_el)
        if count == 5000:
            store_quants(quants, quantdb)
            count = 0
            quants = {}
        formatting.clear_el(consensus_el)
        formatting.clear_el(spectrum)
    return quantdb.fn


def get_quant_data(cons_el):
    """Gets quant data from consensusXML element"""
    quant_out = {}
    for reporter in cons_el.findall('element'):
        quant_out[reporter.attrib['map']] = reporter.attrib['it']
    return quant_out


def store_quants(quants, quantdb):
    sql_quants = []
    for specfile, scannr, quantdata in list(quants.items()):
        for quantmap in sorted(quantdata.keys()):
            sql_quants.append((scannr, quantmap, quantdata[quantmap]))
    quantdb.store_quants(sql_quants)
