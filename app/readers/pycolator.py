from lxml import etree
from . import basereader


def get_root_el(fn):
    rootgen = etree.iterparse(fn, events=('start',))
    return next(rootgen)[1]


def get_namespace(fn):
    root = get_root_el(fn)
    ns = {}
    for prefix in root.nsmap:
        separator = ':'
        if prefix is None:
            separator = None
        ns['xmlns{0}{1}'.format(separator, prefix)] = root.nsmap[prefix]
    return ns


def get_percolator_static_xml(fn, ns):
    #rootgen = etree.iterparse(fn, tag='{%s}percolator_output' % ns['xmlns'],
    #                          events=('start',))
    root = get_root_el(fn)
    #for child in root.getchildren():
    #    root.remove(child)
    process = etree.iterparse(fn, tag='{%s}process_info' % ns['xmlns'],
                              events=('start',))
    root.append(next(process)[1])
    return root


# Stringified element generators interfaces
def generate_psms_multiple_fractions_strings(input_files, ns):
    return basereader.generate_tags_multiple_files_strings(input_files,
                                                           ns, 'psm',
                                                           ['peptide',
                                                            'protein'])


def generate_peptides_multiple_fractions_strings(input_files, ns):
    return basereader.generate_tags_multiple_files_strings(input_files,
                                                           ns, 'peptide',
                                                           ['psm', 'protein'])


# Element generators interfaces
def generate_psms_multiple_fractions(input_files, ns):
    return basereader.generate_tags_multiple_files(input_files, 'psm',
                                                   ['peptide', 'protein'], ns)


def generate_peptides_multiple_fractions(input_files, ns):
    return basereader.generate_tags_multiple_files(input_files, 'peptide',
                                                   ['psm', 'protein'], ns)


def generate_psms(fn, ns):
    return basereader.generate_tags_multiple_files(fn, 'psm',
                                                   ['peptide', 'protein'], ns)


def generate_peptides(fn, ns):
    return basereader.generate_xmltags(fn, 'peptide', ['psm', 'protein'], ns)
