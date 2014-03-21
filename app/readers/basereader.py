from lxml import etree
from app import formatting


def generate_tags_multiple_files(input_files, tag, ignore_tags, ns=None):
    """
    Base generator for percolator xml psm, peptide, protein output.
    """
    if ns is None:
        xmlns = ''
    else:
        xmlns = '{%s}' % ns['xmlns']
    for fn in input_files:
        for ac, el in etree.iterparse(fn):
            if el.tag == '{0}{1}'.format(xmlns, tag):
                yield el
            elif el.tag in ['{0}{1}'.format(xmlns, x) for x in ignore_tags]:
                formatting.clear_el(el)


def generate_tags_multiple_files_strings(input_files, ns, tag, ignore_tags):
    """
    Creates stringified xml output of elements with certain tag.
    ignore_tags are the ones that are not cleared when met
    """
    for el in generate_tags_multiple_files(input_files, tag, ignore_tags, ns):
        yield formatting.string_and_clear(el, ns)


