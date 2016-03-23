import os
from app.readers import xml as basereader


def specfn_quant_generator(specfiles, quantfiles, tag, ignore_tags):
    """Generates tuples of specfile and quant element for general formats"""
    for specfn, qfn in zip(specfiles, quantfiles):
        for quant_el in basereader.generate_xmltags(qfn, tag, ignore_tags):
            yield os.path.basename(specfn), quant_el


def mzmlfn_cons_el_generator(specfiles, consfiles):
    """Returns generation of tuples of spectra file and
    consensusXML quant elements"""
    return specfn_quant_generator(specfiles, consfiles, 'consensusElement',
                                  ['consensusElementList'])


def mzmlfn_feature_generator(specfiles, featfiles):
    """Returns tuple of spectrafile and features of OpenMS
    feature XML format"""
    return specfn_quant_generator(specfiles, featfiles, 'feature',
                                  ['featureList'])


def get_consxml_rt(cons_el):
    """Returns consensusXML in seconds"""
    return cons_el.find('centroid').attrib['rt']


def get_feature_info(feature):
    """Returns a dict with feature information"""
    dimensions = feature.findall('position')
    for dim in dimensions:
        if dim.attrib['dim'] == '0':
            rt = dim.text
        elif dim.attrib['dim'] == '1':
            mz = dim.text
    return {'rt': float(rt), 'mz': float(mz),
            'charge': int(feature.find('charge').text),
            'intensity': float(feature.find('intensity').text),
            }


def get_quantmap(consfile):
    """Returns map of isobaric quant channel names as keys, channel numbers
    values from consensus XML.
    E.g. {'113': '0',
          '114': '1',}
    """
    quantmap = {}
    maplist = basereader.get_element(consfile, 'mapList')
    for mapitem in maplist:
        quantmap[mapitem.attrib['label']] = mapitem.attrib['id']
    return quantmap
