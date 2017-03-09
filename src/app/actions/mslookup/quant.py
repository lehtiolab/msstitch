from decimal import Decimal

from app.readers import openms as openmsreader
DB_STORE_CHUNK = 500000
FEATURE_ALIGN_WINDOW_AMOUNT = 1000
PROTON_MASS = 1.0072


def create_isobaric_quant_lookup(quantdb, specfn_consensus_els, channelmap):
    """Creates an sqlite lookup table of scannrs with quant data.

    spectra - an iterable of tupled (filename, spectra)
    consensus_els - a iterable with consensusElements"""
    # store quantchannels in lookup and generate a db_id vs channel map
    channels_store = ((name,) for name, c_id
                      in sorted(channelmap.items(), key=lambda x: x[1]))
    quantdb.store_channelmap(channels_store)
    channelmap_dbid = {channelmap[ch_name]: ch_id for ch_id, ch_name in
                       quantdb.get_channelmap()}
    quants = []
    mzmlmap = quantdb.get_mzmlfile_map()
    for specfn, consensus_el in specfn_consensus_els:
        rt = openmsreader.get_consxml_rt(consensus_el)
        rt = round(float(Decimal(rt) / 60), 12)
        qdata = get_quant_data(consensus_el)
        spectra_id = quantdb.get_spectra_id(mzmlmap[specfn],
                                            retention_time=rt)
        for channel_no in sorted(qdata.keys()):
            quants.append((spectra_id, channelmap_dbid[channel_no],
                           qdata[channel_no]))
            if len(quants) == DB_STORE_CHUNK:
                quantdb.store_isobaric_quants(quants)
    quantdb.store_isobaric_quants(quants)
    quantdb.index_isobaric_quants()


def create_precursor_quant_lookup(quantdb, mzmlfn_feats, quanttype,
                                  rttol, mztol, mztoltype):
    """Fills quant sqlite with precursor quant from:
        features - generator of xml features from openms
    """
    featparsermap = {'kronik': kronik_featparser,
                     'openms': openms_featparser,
                     }
    features = []
    mzmlmap = quantdb.get_mzmlfile_map()
    for specfn, feat_element in mzmlfn_feats:
        feat = featparsermap[quanttype](feat_element)
        features.append((mzmlmap[specfn], feat['rt'], feat['mz'],
                         feat['charge'], feat['intensity'])
                        )
        if len(features) == DB_STORE_CHUNK:
            quantdb.store_ms1_quants(features)
            features = []
    quantdb.store_ms1_quants(features)
    quantdb.index_precursor_quants()
    align_quants_psms(quantdb, rttol, mztol, mztoltype)


def get_minmax(center, tolerance, toltype=None):
    center = float(center)
    if toltype == 'ppm':
        tolerance = int(tolerance) / 1000000 * center
    else:
        tolerance = float(tolerance)
    return center - tolerance, center + tolerance


def align_quants_psms(quantdb, rt_tolerance, mz_tolerance, mz_toltype):
    allspectra = quantdb.get_spectra_mz_sorted()
    featwindow_max_mz = -1
    spec_feat_store = []
    for spec_id, fn_id, charge, mz, rt in allspectra:
        minmz, maxmz = get_minmax(mz, mz_tolerance, mz_toltype)
        if maxmz > featwindow_max_mz:
            feat_map, featwindow_max_mz = get_precursors_from_window(quantdb,
                                                                     minmz)
        best_feat_id = align_psm(mz, rt, fn_id, charge, feat_map, rt_tolerance)
        if not best_feat_id:
            continue
        spec_feat_store.append((spec_id, best_feat_id))
        if len(spec_feat_store) > DB_STORE_CHUNK:
            quantdb.store_ms1_alignments(spec_feat_store)
            spec_feat_store = []
    quantdb.store_ms1_alignments(spec_feat_store)


def align_psm(psm_mz, psm_rt, fn_id, charge, featmap, rttol):
    minrt, maxrt = get_minmax(psm_rt, rttol / 60)
    alignments = {}
    try:
        featlist = featmap[fn_id][charge]
    except KeyError:
        return False
    for feat_mz, feat_rt, feat_id in featlist:
        if abs(psm_rt - feat_rt) > rttol:
            continue
        alignments[abs(psm_mz - feat_mz)] = feat_id
    try:
        return alignments[min(alignments)]
    except ValueError:
        return False


def get_precursors_from_window(quantdb, minmz):
    """Returns a dict of a specified amount of features from the
    ms1 quant database, and the highest mz of those features"""
    featmap = {}
    mz = False
    features = quantdb.get_precursor_quant_window(FEATURE_ALIGN_WINDOW_AMOUNT,
                                                  minmz)
    for feat_id, fn_id, charge, mz, rt in features:
        try:
            featmap[fn_id][charge].append((mz, rt, feat_id))
        except KeyError:
            try:
                featmap[fn_id][charge] = [(mz, rt, feat_id)]
            except KeyError:
                featmap[fn_id] = {charge: [(mz, rt, feat_id)]}
    return featmap, mz


def kronik_featparser(feature):
    charge = int(feature['Charge'])
    mz = (float(feature['Monoisotopic Mass']) + charge * PROTON_MASS) / charge
    return {'rt': round(float(feature['Best RTime']), 12),
            'mz': mz,
            'charge': charge,
            'intensity': float(feature['Best Intensity']),
            }


def openms_featparser(feature):
    feat = openmsreader.get_feature_info(feature)
    feat['rt'] = round(float(Decimal(feat['rt']) / 60), 12)
    return feat


def get_quant_data(cons_el):
    """Gets quant data from consensusXML element"""
    quant_out = {}
    for reporter in cons_el.findall('.//element'):
        quant_out[reporter.attrib['map']] = reporter.attrib['it']
    return quant_out
