import sys
from lxml import etree
import readers

def stringify_strip_namespace_declaration(el, ns):
    strxml = etree.tostring(el)
    strxml = strxml.replace('xmlns="{0}" '.format(ns['xmlns']), '')
    strxml = strxml.replace('xmlns:p="{0}" '.format(ns['xmlns:p']), '')
    strxml = strxml.replace('xmlns:xsi="{0}" '.format(ns['xmlns:xsi']), '')
    return strxml

def target_decoy_generator(element_generator, decoy, ns):
    for ev,el in element_generator:
        if el.attrib['{%s}decoy' % ns['xmlns']] == decoy:
            strxml = stringify_strip_namespace_declaration(el, ns)
        else:
            continue
        yield strxml


def split_target_decoy(fn, ns):
    split_elements = {'target': {}, 'decoy': {}}
    feats_to_process = ['psm', 'peptide']
    for feat in feats_to_process:
        targetfeatgen = etree.iterparse(fn, tag='{%s}%s' % (ns['xmlns'], feat))
        decoyfeatgen = etree.iterparse(fn, tag='{%s}%s' % (ns['xmlns'], feat))
        split_elements['target'][feat] = target_decoy_generator(targetfeatgen, 'false', ns)
        split_elements['decoy'][feat] = target_decoy_generator(decoyfeatgen, 'true', ns)
        
    return split_elements

def filter_unique_peptides(input_files, score, ns):
    """ Filters unique peptides from multiple Percolator output XML files.
        Takes a dir with a set of XML files, a score to filter on and a namespace.
        Outputs an ElementTree.
    """
    if len(input_files) == 1:
        print >> sys.stdout, 'WARNING: You are attempting to filter unique peptides from a non-fractioned dataset. This may be unneccessary.'
    scores = {'q':'q_value', 'pep':'pep', 'p':'p_value', 'svm':'svm_score'}
    pepgen = \
        readers.generate_peptides_multiple_fractions(input_files, ns)
    highest = {}
    for el in pepgen:
        featscore = float(el.xpath('xmlns:%s' % scores[score], namespaces=ns)[0].text)
        seq = el.attrib['{%s}peptide_id' % ns['xmlns']]

        if seq not in highest:
            highest[seq] = {'pep_el': el, 'score': featscore}
        if score == 'svm': # greater than score is accepted
            if featscore > highest[seq]['score']:
                highest[seq] = {'pep_el': el, 'score': featscore}
        else: # lower than score is accepted
            if featscore < highest[seq]['score']:
                highest[seq] = {'pep_el': el, 'score': featscore}
    
    for pep in highest.values():
        yield stringify_strip_namespace_declaration(pep['pep_el'], ns)
