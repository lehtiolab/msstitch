from lxml import etree
import formatting

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
    return generate_tags_multiple_files_strings(input_files, ns, 'psm',
                                                ['peptide', 'protein'])

def generate_peptides_multiple_fractions_strings(input_files, ns):
    return generate_tags_multiple_files_strings(input_files, ns, 'peptide',
                                                ['psm', 'protein'])

# Element generators interfaces
def generate_psms_multiple_fractions(input_files, ns):
    return generate_tags_multiple_files(input_files, ns, 'psm',
                                                ['peptide', 'protein'])

def generate_peptides_multiple_fractions(input_files, ns):
    return generate_tags_multiple_files(input_files, ns, 'peptide',
                                                ['psm', 'protein'])


# Generators that actually do the work
def generate_tags_multiple_files_strings(input_files, ns, tag, ignore_tags):
    """
    Creates stringified output for percolator elements of certain tag.
    """
    for el in generate_tags_multiple_files(input_files, ns, tag, ignore_tags):
        yield formatting.string_and_clear(el, ns)

def generate_tags_multiple_files(input_files, ns, tag, ignore_tags):
    """
    Base generator for percolator xml psm, peptide, protein output.
    """
    for fn in input_files:
        for ac,el in etree.iterparse(fn):
            if el.tag=='{%s}%s' % (ns['xmlns'], tag):
                yield el
            elif el.tag in ['{%s}%s' % (ns['xmlns'], x) for x in ignore_tags]:
                formatting.clear_el(el)



def generate_peptides_by_seq_multiple_fractions(input_files, seq, ns):
    # Obsolete, deprecate? FIXME
    for fn in input_files:
        for ac,el in etree.iterparse(fn, tag='{%s}peptide' % ns['xmlns']):
            if el.attrib['{%s}peptide_id' % ns['xmlns']] != seq:
                continue
            yield el


