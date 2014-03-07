import basereader


def psm_generator(mzidfile):
    return basereader.generate_tags_multiple_files(
        [mzidfile],
        'SpectrumIdentificationResult',
        ['MzIdentML',
         'DataCollection',
         'AnalysisData',
         'SpectrumIdentificationList'])


def mzml_generator(mzmlfile):
    return basereader.generate_tags_multiple_files(
        [mzmlfile],
        'spectrum',
        ['indexedmzML',
         'mzML',
         'run',
         'spectrumList'])


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


def get_scan_nr(element, attribname):
    """General method to get a scan nr from xml element of mzML/mzIdentML"""
    info = element.attrib[attribname]
    infomap = {y[0]: y[1] for x in info.split() for y in x.split('=')}
    return infomap['scan']


def get_spec_scan_nr(spectrum):
    """Returns scan number of mzML spectrum as a str."""
    return get_scan_nr(spectrum, 'id')


def get_psm_scan_nr(psm):
    """Returns scan nr of an mzIdentML PSM as a str. The PSM is given
    as a SpectrumIdentificationResult element."""
    return get_scan_nr(psm, 'spectrumID')
