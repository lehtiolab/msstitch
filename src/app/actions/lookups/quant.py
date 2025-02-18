import os
from decimal import Decimal
from collections import defaultdict

from app.readers import openms as openmsreader
DB_STORE_CHUNK = 500000
FEATURE_ALIGN_WINDOW_AMOUNT = 1000
PROTON_MASS = 1.0072


def create_isobaric_quant_lookup(quantdb, specfns, isobaricfns):
    """Creates an sqlite lookup table of scannrs with quant data.

    spectra - an iterable of tupled (filename, spectra)
    consensus_els - a iterable with consensusElements"""
    mzmlmap = quantdb.get_mzmlfile_map()
    for specfpath, isofn in zip(specfns, isobaricfns):
        specfn = os.path.basename(specfpath)
        quants, pifs = [], []
        # store quantchannels in lookup and generate a db_id vs channel map
        channelmap = openmsreader.get_quantmap(isofn)
        channels_store = ((name,) for name, c_id
                          in sorted(channelmap.items(), key=lambda x: int(x[1])))
        quantdb.store_channelmap(channels_store)
        channelmap_dbid = {channelmap[ch_name]: ch_id for ch_id, ch_name in
                quantdb.get_channelmap() if ch_name in channelmap}
        # Get DB spectra IDs keyed by scan SIDs
        specmap = quantdb.get_specmap(mzmlmap[specfn])
        for _, consensus_el in openmsreader.mzmlfn_cons_el_generator([specfn], [isofn]):
            pif, sid = openmsreader.get_consxml_sid_pif(consensus_el)
            # FIXME what if sid not in specmap? error then
            spectra_id = specmap[sid]['id']
            if pif:
                pifs.append((spectra_id, float(pif)))
            qdata = get_quant_data(consensus_el)
            for channel_no in sorted(qdata.keys()):
                quants.append((spectra_id, channelmap_dbid[channel_no],
                               qdata[channel_no]))
                if len(quants) == DB_STORE_CHUNK:
                    quantdb.store_isobaric_quants(quants, pifs)
                    quants, pifs = [], []
        quantdb.store_isobaric_quants(quants, pifs)
    quantdb.index_isobaric_quants()


def create_precursor_quant_lookup(quantdb, mzmlfn_feats, sum_or_apex, quanttype,
        rttol, mztol, mztoltype):
    """Fills quant sqlite with precursor quant from:
        features - generator of xml features from openms
    """
    featparsermap = {'kronik': kronik_featparser,
                     'openms': openms_featparser,
                     'dinosaur': dinosaur_featparser,
                     }
    features, fwhms = [], []
    mzmlmap = quantdb.get_mzmlfile_map()
    for specfn, feat_element in mzmlfn_feats:
        feat = featparsermap[quanttype](feat_element, sum_or_apex)
        features.append((mzmlmap[specfn], feat['rt'], feat['mz'],
                         feat['charge'], feat['intensity'])
                        )
        fwhms.append(feat['fwhm'])
        if len(features) == DB_STORE_CHUNK:
            id_feats = quantdb.store_ms1_quants(features)
            if quanttype == 'dinosaur':
                quantdb.store_fwhm(zip([x[0] for x in id_feats], fwhms))
            features, fwhms = [], []
    id_feats = quantdb.store_ms1_quants(features)
    if quanttype == 'dinosaur':
        quantdb.store_fwhm(zip([x[0] for x in id_feats], fwhms))
    quantdb.index_precursor_quants()
    align_quants_psms(quantdb, rttol, mztol, mztoltype)
    quantdb.index_aligned_quants()


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
    active_fn = None
    for spec_id, fn_id, charge, mz, rt in allspectra:
        if fn_id != active_fn:
            active_fn = fn_id
            featwindow_max_mz = -1
            fnfeats = quantdb.get_fnfeats(fn_id)
        minmz, maxmz = get_minmax(mz, mz_tolerance, mz_toltype)
        if maxmz > featwindow_max_mz:
            feat_map, featwindow_max_mz = get_precursors_from_window(fnfeats,
                                                                     minmz)
        best_feat_id = align_psm(mz, rt, charge, feat_map, rt_tolerance)
        if not best_feat_id:
            continue
        spec_feat_store.append((spec_id, best_feat_id))
        if len(spec_feat_store) > DB_STORE_CHUNK:
            quantdb.store_ms1_alignments(spec_feat_store)
            spec_feat_store = []
    quantdb.store_ms1_alignments(spec_feat_store)


def align_psm(psm_mz, psm_rt, charge, featmap, rttol):
    rttol_min = rttol / 60
    alignments = defaultdict(list)
    try:
        featlist = featmap[charge]
    except KeyError:
        return False
    for feat_mz, feat_rt, feat_id in featlist:
        rtdiff = abs(psm_rt - feat_rt)
        mzdiff = abs(psm_mz - feat_mz)
        if rtdiff > rttol_min:
            continue
        alignments[mzdiff].append((feat_id, rtdiff))
    try:
        best_mz = alignments[min(alignments)]
    except ValueError:
        return False
    if len(best_mz) > 1:
        print('Multiple features of same m/z found within RT window, possibly '
                'your RT tolerance is set too large')
        return sorted(best_mz, key=lambda x: x[1])[0][0]
    else:
        return best_mz[0][0]


def get_precursors_from_window(mzfeatmap, minmz):
    """Returns a dict of a specified amount of features from the
    ms1 quant database, and the highest mz of those features"""
    chargemap = {}
    mz, amount_feats = False, 0
    features = []
    for mz, feat_id, charge, rt in mzfeatmap:
        if mz < minmz:
            continue
        elif amount_feats > FEATURE_ALIGN_WINDOW_AMOUNT:
            break
        amount_feats += 1
        try:
            chargemap[charge].append((mz, rt, feat_id))
        except KeyError:
            chargemap[charge] = [(mz, rt, feat_id)]
    return chargemap, mz


def kronik_featparser(feature, sum_or_apex):
    charge = int(feature['Charge'])
    mz = (float(feature['Monoisotopic Mass']) + charge * PROTON_MASS) / charge
    intkey = {'sum': 'Summed Intensity', 'apex': 'Best Intensity'}[sum_or_apex]
    return {'rt': round(float(feature['Best RTime']), 12),
            'mz': mz,
            'charge': charge,
            'intensity': float(feature[intkey]),
            'fwhm': False,
            }


def dinosaur_featparser(feature, sum_or_apex):
    # FIXME is kronik, change, add FWHM
    #mz = (float(feature['mz']) + charge * PROTON_MASS) / charge
    intkey = {'sum': 'intensitySum', 'apex': 'intensityApex'}[sum_or_apex]
    return {'rt': round(float(feature['rtApex']), 12),
            'mz': float(feature['mz']),
            'charge': int(feature['charge']),
            'intensity': float(feature[intkey]),
            'fwhm': float(feature['fwhm']),
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
