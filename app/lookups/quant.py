from app import sqlite
from app import formatting
from app.readers import openms as readers


def create_quant_lookup(fn_spectra, consensus_els):
    """Creates an sqlite lookup table of scannrs with quant data.

    spectra - an iterable of tupled (filename, spectra)
    consensus_els - a iterable with consensusElements"""
    quantdb = sqlite.QuantDB()
    quantdb.create_quantdb()
    quants = {}
    count = 0
    for consensus_el in consensus_els:
        count += 1
        rt = readers.get_consxml_rt(consensus_el)
        fn, spec_scan_nr = get_spec_scan_nr(fn_spectra, rt)
        qdata = get_quant_data(consensus_el)
        try:
            quants[fn][spec_scan_nr] = qdata
        except KeyError:
            quants[fn] = {}
            quants[fn][spec_scan_nr] = qdata
        if count == 5000:
            store_quants(quants, quantdb)
            count = 0
            quants = {}
        formatting.clear_el(consensus_el)
    store_quants(quants, quantdb)
    return quantdb.fn


def get_spec_scan_nr(fn_spectra, cons_rt):
    """Returns mzML filename and scan nr of a spectrum that belongs to the
    retention time in cons_rt, when given a generator of multiple
    spectra files in fn_spectra.
    The generator is a tuple of fn, spectrum, namespace."""
    rt = round(cons_rt, 8)
    for fn, spectrum, ns in fn_spectra:
        if rt == round(readers.get_mzml_rt(spectrum, ns), 8):
            scan_nr = readers.get_spec_scan_nr(spectrum)
            return fn, scan_nr
        formatting.clear_el(spectrum)
    return False


def get_quant_data(cons_el):
    """Gets quant data from consensusXML element"""
    quant_out = {}
    for reporter in cons_el.findall('.//element'):
        quant_out[reporter.attrib['map']] = reporter.attrib['it']
    return quant_out


def store_quants(quants, quantdb):
    sql_quants = []
    for specfn in quants:
        for scannr, quantdata in quants[specfn].items():
            for quantmap in sorted(quantdata.keys()):
                sql_quants.append((specfn, scannr,
                                   quantmap, quantdata[quantmap]))
    quantdb.store_quants(sql_quants)
