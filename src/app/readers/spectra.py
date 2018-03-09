import os
from app.readers import xml as basereader
from app.readers import ml
from app.readers import xmlformatting as formatting


def mzmlfn_ms2_spectra_generator(mzmlfiles):
    for fn, spec, ns in mzmlfn_spectra_generator(mzmlfiles):
        specparams = get_all_cvparams(spec, ns)
        mslvl = fetch_cvparam_value_by_name(specparams, 'ms level')
        if mslvl != '2':
            continue
        scannr = get_spec_scan_nr(spec)
        rt = fetch_cvparams_values_from_subel(spec, 'scan',
                                              ['scan start time'], ns)
        iit = fetch_cvparams_values_from_subel(spec, 'scan',
                                               ['ion injection time'], ns)
        mz, charge = fetch_cvparams_values_from_subel(spec, 'selectedIon',
                                                      ['selected ion m/z',
                                                       'charge state'], ns)
        yield fn, {'scan': scannr, 'rt': rt[0], 'iit': iit[0], 'mz': mz,
                   'charge': charge}
        formatting.clear_el(spec)


def mzmlfn_spectra_generator(mzmlfiles):
    for fn in mzmlfiles:
        ns = basereader.get_namespace(fn)
        spectra = basereader.generate_xmltags(
            fn,
            'spectrum',
            ['offset',
             ],
            ns)
        for spectrum in spectra:
            yield os.path.basename(fn), spectrum, ns


def fetch_cvparams_values_from_subel(base, subelname, paramnames, ns):
    """Searches a base element for subelement by name, then takes the
    cvParams of that subelement and returns the values as a list
    for the paramnames that match. Value order in list equals input
    paramnames order."""
    sub_el = basereader.find_element_xpath(base, subelname, ns)
    cvparams = get_all_cvparams(sub_el, ns)
    output = []
    for param in paramnames:
        output.append(fetch_cvparam_value_by_name(cvparams, param))
    return output


def fetch_cvparam_value_by_name(params, name):
    for cvparam in params:
        try:
            if cvparam.attrib['name'] == name:
                return cvparam.attrib['value']
        except KeyError:
            pass
    return False


def get_all_cvparams(element, ns):
    return element.findall('{%s}cvParam' % ns['xmlns'])


def get_spec_scan_nr(spectrum):
    """Returns scan number of mzML spectrum as a str."""
    return ml.get_scan_nr(spectrum, 'id')
