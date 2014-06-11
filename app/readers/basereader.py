from lxml import etree
import itertools
from app import formatting


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


def generate_xmltags(fn, tag, ignore_tags, ns=None):
    """
    Base generator for percolator xml psm, peptide, protein output,
    as well as for mzML, mzIdentML.
    ignore_tags are the ones that are not cleared when met by parser.
    """
    if ns is None:
        xmlns = ''
    else:
        xmlns = '{%s}' % ns['xmlns']
    for ac, el in etree.iterparse(fn):
        if el.tag == '{0}{1}'.format(xmlns, tag):
            yield el
        elif el.tag in ['{0}{1}'.format(xmlns, x) for x in
                        ignore_tags]:
            formatting.clear_el(el)


def generate_tsv_lines_multifile(fns):
    return itertools.chain.from_iterable([generate_tsv_psms(fn)
                                          for fn in fns])


def generate_tsv_psms(fn):
    """Returns dicts with header-keys and psm statistic values"""
    with open(fn) as fp:
        header = next(fp)
        for line in fp:
            yield {x: line.strip().split('\t')[x] for x in header}
