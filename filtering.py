import re
from lxml import etree

def get_peptide_seq(peptide, ns):
    return peptide.attrib['{%s}peptide_id' % ns['xmlns']]

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


def filter_known_searchspace(peptides, searchspace, ns):
    """Yields peptides from generator as long as their sequence is not found in
    known search space dict. Useful for excluding peptides that are found in
    e.g. ENSEMBL or similar"""
    for peptide in peptides:
        seq = get_peptide_seq(peptide, ns)
        # Loose modifications
        seq = re.sub('\[UNIMOD:\d*\]', '', seq)
        try:
            searchspace[seq]
        except KeyError:
            yield peptide
        else:
            pass

def filter_unique_peptides(peptides, score, ns):
    """ Filters unique peptides from multiple Percolator output XML files.
        Takes a dir with a set of XML files, a score to filter on and a namespace.
        Outputs an ElementTree.
    """
    scores = {'q':'q_value', 'pep':'pep', 'p':'p_value', 'svm':'svm_score'}
    highest = {}
    for el in peptides:
        featscore = float(el.xpath('xmlns:%s' % scores[score], namespaces=ns)[0].text)
        seq = get_peptide_seq(el, ns)
        
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
