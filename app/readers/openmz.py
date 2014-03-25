from . import basereader
from . import ml


def mzml_generator(mzmlfiles):
    for fn in mzmlfiles:
        spectra = basereader.generate_tags_multiple_files(
            fn,
            'spectrum',
            ['indexedmzML',
             'mzML',
             'run',
             'spectrumList'])
        for spectrum in spectra:
            yield fn, spectrum


def quant_generator(consfile):
    return basereader.generate_tags_multiple_files(
        [consfile],
        'consensusElement',
        ['consensusXML', 'consensusElementList'])


def get_consxml_rt(cons_el):
    """Returns consensusXML in minutes"""
    rt = cons_el.find('centroid').attrib['rt']
    return int(rt) / 60


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
