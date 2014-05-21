from . import basereader
from . import ml


def mzml_generator(mzmlfiles):
    for fn in mzmlfiles:
        ns = basereader.get_namespace(fn)
        spectra = basereader.generate_xmltags(
            fn,
            'spectrum',
            ['run',
             'spectrumList'],
            ns)
        for spectrum in spectra:
            yield fn, spectrum, ns


def quant_generator(consfile):
    return basereader.generate_xmltags(
        consfile,
        'consensusElement',
        ['consensusElementList'],
        )


def get_consxml_rt(cons_el):
    """Returns consensusXML in minutes"""
    rt = cons_el.find('centroid').attrib['rt']
    return float(rt) / 60


def get_mzml_rt(spectrum, ns):
    """Return RT from mzml spectrum element, in minutes"""
    scan = spectrum.find('.//{%s}scan' % ns['xmlns'])
    for cvparam in scan.findall('{%s}cvParam' % ns['xmlns']):
        try:
            if cvparam.attrib['name'] == 'scan start time':
                return float(cvparam.attrib['value'])
        except KeyError:
            pass


def get_spec_scan_nr(spectrum):
    """Returns scan number of mzML spectrum as a str."""
    return ml.get_scan_nr(spectrum, 'id')
