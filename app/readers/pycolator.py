from lxml import etree
from . import basereader


def get_namespace(fn):
    ns = {'xmlns': 'http://per-colator.com/percolator_out/14',
          'xmlns:p': 'http://per-colator.com/percolator_out/14',
          'xmlns:xsi': 'http://www.w3.org/2001/XMLSchema-instance',
          }
    return ns
    # FIXME lookup from file


def get_percolator_static_xml(fn, ns):
    rootgen = etree.iterparse(fn, tag='{%s}percolator_output' % ns['xmlns'],
                              events=('start',))
    root = next(rootgen)[1]
    for child in root.getchildren():
        root.remove(child)
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
