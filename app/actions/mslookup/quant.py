from decimal import Decimal, getcontext

from app.readers import openms as openmsreader


def create_isobaric_quant_lookup(quantdb, specfn_consensus_els, channelmap):
    """Creates an sqlite lookup table of scannrs with quant data.

    spectra - an iterable of tupled (filename, spectra)
    consensus_els - a iterable with consensusElements"""
    quants = []
    mzmlmap = quantdb.get_mzmlfile_map()
    for specfn, consensus_el in specfn_consensus_els:
        rt = openmsreader.get_consxml_rt(consensus_el)
        rt = float(Decimal(rt) / 60)
        qdata = get_quant_data(consensus_el)
        spectra_id = quantdb.get_spectra_id(mzmlmap[specfn], rt)
        for channel_no in sorted(qdata.keys()):
            quants.append((spectra_id, channelmap[channel_no],
                           qdata[channel_no]))
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
    mzmlmap = quantdb.get_mzmlfile_map()
    for specfn, feat_element in mzmlfn_featsxml:
        feat = openmsreader.get_feature_info(feat_element)
        feat['rt'] = float(Decimal(feat['rt']) / 60)
        features.append((mzmlmap[specfn], feat['rt'], feat['mz'],
                         feat['charge'], feat['intensity'])
                        )
        if len(features) == 5000:
            quantdb.store_ms1_quants(features)
            features = []
    quantdb.store_ms1_quants(features)
    quantdb.index_precursor_quants()


def get_quant_data(cons_el):
    """Gets quant data from consensusXML element"""
    quant_out = {}
    for reporter in cons_el.findall('.//element'):
        quant_out[reporter.attrib['map']] = reporter.attrib['it']
    return quant_out
