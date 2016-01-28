from app.readers import xmlformatting as formatting
from app.readers import pycolator as readers


def get_score(elements, ns, scoretype='svm_score'):
    for el in elements:
        score = el.xpath('xmlns:{0}'.format(scoretype), namespaces=ns)[0].text
        formatting.clear_el(el)
        yield score


def prepare_qvality_input(fn, feattype, get_static_percolator_output):
    featextractors = {'peptide': readers.generate_peptides,
                      'psm': readers.generate_psms
                      }

    ns, static_xml = get_static_percolator_output(fn)
    features_for_qvality = featextractors[feattype](fn, ns)
    return get_score(features_for_qvality, ns)
