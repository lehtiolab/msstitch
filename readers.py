import sys
from lxml import etree
import filtering

# Deprecate? FIXME
def readPercout(fname):
    doc = None
    try:
        doc = etree.parse(fname)
    except Exception:
        sys.stderr.write('Could not parse XML provided in %s or error reading file. \n' % (fname))
    return doc

def get_namespace(fn):
    ns = {'xmlns':'http://per-colator.com/percolator_out/14',
    'xmlns:p' : 'http://per-colator.com/percolator_out/14',
    'xmlns:xsi':'http://www.w3.org/2001/XMLSchema-instance',
    }
    return ns
    # FIXME lookup from file


def get_percolator_static_xml(fn, ns):
    rootgen = etree.iterparse(fn, tag='{%s}percolator_output' % ns['xmlns'], events=('start',))
    root = rootgen.next()[1]
    for child in root.getchildren():
        root.remove(child)
    process = etree.iterparse(fn, tag='{%s}process_info' % ns['xmlns'], events=('start',))
    root.append(process.next()[1])
    return root



# Stringified element generators interfaces
def generate_psms_multiple_fractions_strings(input_files, ns):
    return generate_tags_multiple_files_strings(input_files, ns, 'psm')

def generate_peptides_multiple_fractions_strings(input_files, ns):
    return generate_tags_multiple_files_strings(input_files, ns, 'peptide')

# Element generators
def generate_psms_multiple_fractions(input_files, ns):
    return generate_tags_multiple_files(input_files, ns, 'psm')

def generate_peptides_multiple_fractions(input_files, ns):
    return generate_tags_multiple_files(input_files, ns, 'peptide')


# String and element generators
def generate_tags_multiple_files_strings(input_files, ns, tag):
    """
    Creates stringified output for percolator elements of certain tag.
    """
    for el in generate_tags_multiple_files(input_files, ns, tag):
        str_el = filtering.stringify_strip_namespace_declaration(el, ns)
        filtering.clear_el(el)
        yield str_el
        
def generate_tags_multiple_files(input_files, ns, tag):
    """
    Base generator for percolator xml psm, peptide, protein output.
    """
    for fn in input_files:
        for ac,el in etree.iterparse(fn):
            if el.tag=='{%s}%s' % (ns['xmlns'], tag):
                yield el
            else:
                filtering.clear_el(el)



def generate_peptides_by_seq_multiple_fractions(input_files, seq, ns):
    # Obsolete, deprecate? FIXME
    for fn in input_files:
        for ac,el in etree.iterparse(fn, tag='{%s}peptide' % ns['xmlns']):
            if el.attrib['{%s}peptide_id' % ns['xmlns']] != seq:
                continue
            yield el


