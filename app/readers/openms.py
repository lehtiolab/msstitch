import os
import app.readers.xml as basereader


def quant_generator(consfiles):
    return basereader.generate_tags_multiple_files(
        consfiles,
        'consensusElement',
        ['consensusElementList'],
    )


def mzmlfn_feature_generator(spectrafiles, featfiles):
    """Returns tuple of spectrafile and features of OpenMS
    feature XML format"""
    for specfn, featfn in zip(spectrafiles, featfiles):
        ns = basereader.get_namespace(featfn)
        features = basereader.generate_xmltags(
            featfn,
            'feature',
            ['featureList'],
            ns)
        for feature in features:
            yield os.path.basename(specfn), feature


def get_consxml_rt(cons_el):
    """Returns consensusXML in minutes"""
    rt = cons_el.find('centroid').attrib['rt']
    return float(rt)


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
