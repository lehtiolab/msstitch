from lxml import etree
import formatting


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


