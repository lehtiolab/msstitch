import basereader


def psm_generator(mzidfile):
    return basereaders.generate_tags_multiple_files(
        [mzidfile],
        'SpectrumIdentificationResult',
        ['MzIdentML',
         'DataCollection',
         'AnalysisData',
         'SpectrumIdentificationList'])


def mzml_generator(mzmlfile):
    return basereaders.generate_tags_multiple_files(
        [mzmlfile],
        'spectrum',
        ['indexedmzML',
         'mzML',
         'run',
         'spectrumList'])


def quant_generator(consfile):
    return basereaders.generate_tags_multiple_files(
        [consfile],
        'consensusElement',
        ['consensusXML', 'consensusElementList'])
