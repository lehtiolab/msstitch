from decimal import Decimal, getcontext

from app.lookups.sqlite import quant as sqlite
from app.readers import openms as openmsreader
from app.readers import spectra as specreader


def initiate_quant_lookup(workdir):
    """Creates a quant db sqlite file and fills a table with spectra mzml
    info on scan nrs and retention times"""
    quantdb = sqlite.QuantDB()
    quantdb.create_quantdb(workdir)
    return quantdb


def get_quant_lookup(quantfn):
    return sqlite.QuantDB(quantfn)


def create_isobaric_quant_lookup(quantdb, specfn_consensus_els):
    """Creates an sqlite lookup table of scannrs with quant data.

    spectra - an iterable of tupled (filename, spectra)
    consensus_els - a iterable with consensusElements"""
    quants = []
    for specfn, consensus_el in specfn_consensus_els:
        rt = openmsreader.get_consxml_rt(consensus_el)
        qdata = get_quant_data(consensus_el)
        for quantmap in sorted(qdata.keys()):
            quants.append((specfn, rt, quantmap, qdata[quantmap]))
            if len(quants) == 5000:
                quantdb.store_isobaric_quants(quants)
    quantdb.store_isobaric_quants(quants)
    quantdb.index_isobaric_quants()


def create_precursor_quant_lookup(quantdb, mzmlfn_featsxml):
    """Fills quant sqlite with precursor quant from:
        features - generator of xml features from openms
    """
    features = []
    getcontext().prec = 14  # sets decimal point precision
    for specfn, feat_element in mzmlfn_featsxml:
        feat = openmsreader.get_feature_info(feat_element)
        feat['rt'] = float(Decimal(feat['rt']))
        features.append((specfn, feat['rt'], feat['mz'],
                         feat['charge'], feat['intensity'])
                        )
        if len(features) == 5000:
            quantdb.store_ms1_quants(features)
            features = []
    quantdb.store_ms1_quants(features)
    quantdb.index_precursor_quants()


def create_spectra_lookup(quantdb, fn_spectra):
    """Stores all spectra rt and scan nr in db"""
    to_store = []
    for fn, spectrum, ns in fn_spectra:
        mzml_rt = float(Decimal(specreader.get_mzml_rt(spectrum, ns)))
        scan_nr = specreader.get_spec_scan_nr(spectrum)
        to_store.append((fn, scan_nr, mzml_rt))
        if len(to_store) == 5000:
            quantdb.store_mzmls(to_store)
            to_store = []
    quantdb.store_mzmls(to_store)


def get_quant_data(cons_el):
    """Gets quant data from consensusXML element"""
    quant_out = {}
    for reporter in cons_el.findall('.//element'):
        quant_out[reporter.attrib['map']] = reporter.attrib['it']
    return quant_out
