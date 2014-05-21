from . import basereader
from . import ml


def mzml_generator(mzmlfiles):
    for fn in mzmlfiles:
        ns = basereader.get_namespace(fn)
        spectra = basereader.generate_xmltags(
            fn,
            'spectrum',
            ['indexedmzML',
             'mzML',
             'run',
             'spectrumList'],
            ns)
        for spectrum in spectra:
            yield fn, spectrum


def quant_generator(consfile, ns):
    return basereader.generate_xmltags(
        [consfile],
        'consensusElement',
        ['consensusXML', 'consensusElementList'],
        ns)


def get_consxml_rt(cons_el):
    """Returns consensusXML in minutes"""
    rt = cons_el.find('centroid').attrib['rt']
    return float(rt) / 60


def get_mzml_rt(spectrum):
    """Return RT from mzml spectrum element, in minutes"""
    scan = spectrum.find('.//scan')
    for cvparam in scan.findall('cvParam'):
        try:
            if cvparam.attrib['name'] == 'scan start time':
                return int(cvparam['value'])
        except KeyError:
            pass


def get_spec_scan_nr(spectrum):
    """Returns scan number of mzML spectrum as a str."""
    return ml.get_scan_nr(spectrum, 'id')
