import sys
from lxml import etree

def readPercout(fname):
    doc = None
    try:
        doc = etree.parse(fname)
    except Exception:
        sys.stderr.write('Could not parse XML provided in %s or error reading file. \n' % (fname))
    return doc

def get_percolator_static_xml(fn, ns):
    rootgen = etree.iterparse(fn, tag='{%s}percolator_output' % ns['xmlns'], events=('start',))
    root = rootgen.next()[1]
    for child in root.getchildren():
        root.remove(child)
    process = etree.iterparse(fn, tag='{%s}process_info' % ns['xmlns'], events=('start',))
    root.append(process.next()[1])
    return root


def get_psms(fn, ns):
    return etree.iterparse(fn, tag='{%s}psm' % ns['xmlns'])
    
def get_peptides_multiple_fractions(input_files, ns):
    peptide_generators_fractions = []
    for fn in input_files:
        peptide_generators_fractions.append(etree.iterparse(fn, tag='{%s}peptide' % ns['xmlns']))
