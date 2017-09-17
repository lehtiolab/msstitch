from lxml import etree
import itertools
from app.readers import xmlformatting as formatting


def get_namespace_from_top(fn, key='xmlns'):
    ac, el = next(etree.iterparse(fn))
    return {'xmlns': el.nsmap[key]}


def get_root_el(fn):
    rootgen = etree.iterparse(fn, events=('start',))
    root = next(rootgen)[1]
    for child in root.getchildren():
        root.remove(child)
    return root


def get_namespace(fn):
    root = get_root_el(fn)
    ns = {}
    for prefix in root.nsmap:
        separator = ':'
        nsprefix = prefix
        if prefix is None:
            nsprefix = ''
            separator = ''
        ns['xmlns{0}{1}'.format(separator, nsprefix)] = root.nsmap[prefix]
    return ns


def find_element_xpath(base, name, ns):
    return base.find('.//{%s}%s' % (ns['xmlns'], name))


def generate_tags_multiple_files(input_files, tag, ignore_tags, ns=None):
    """
    Calls xmltag generator for multiple files.
    """
    return itertools.chain.from_iterable([generate_xmltags(
        fn, tag, ignore_tags, ns) for fn in input_files])


def generate_tags_multiple_files_strings(input_files, ns, tag, ignore_tags):
    """
    Creates stringified xml output of elements with certain tag.
    """
    for el in generate_tags_multiple_files(input_files, tag, ignore_tags, ns):
        yield formatting.string_and_clear(el, ns)


def generate_xmltags(fn, returntag, ignore_tags, ns=None):
    """
    Base generator for percolator xml psm, peptide, protein output,
    as well as for mzML, mzIdentML.
    ignore_tags are the ones that are cleared when met by parser.
    """
    xmlns = create_namespace(ns)
    ns_ignore = ['{0}{1}'.format(xmlns, x) for x in ignore_tags]
    for ac, el in etree.iterparse(fn):
        if el.tag == '{0}{1}'.format(xmlns, returntag):
            yield el
        elif el.tag in ns_ignore:
            formatting.clear_el(el)


def get_element(fn, tag, ns=None):
    xmlns = create_namespace(ns)
    for ac, el in etree.iterparse(fn):
        if el.tag == '{0}{1}'.format(xmlns, tag):
            return el


def create_namespace(ns):
    if ns is None:
        return ''
    else:
        return '{%s}' % ns['xmlns']
