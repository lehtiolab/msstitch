from decimal import Decimal, getcontext

from app import sqlite
from app import formatting
from app.readers import openms as openmsreader
from app.readers import spectra as specreader


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
        rt = openmsreader.get_consxml_rt(consensus_el)
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
    quantdb.index_quants()
    return quantdb.fn


def get_spec_scan_nr(fn_spectra, cons_rt):
    """Returns mzML filename and scan nr of a spectrum that belongs to the
    retention time in cons_rt, when given a generator of multiple
    spectra files in fn_spectra.
    The generator is a tuple of fn, spectrum, namespace."""
    getcontext().prec = 8  # sets decimal point precision
    cons_rt_dec = Decimal(cons_rt)
    for fn, spectrum, ns in fn_spectra:
        mzml_rt = Decimal(specreader.get_mzml_rt(spectrum, ns))
        if float(cons_rt_dec/60) == float(mzml_rt):
            scan_nr = specreader.get_spec_scan_nr(spectrum)
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
