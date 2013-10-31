import sys
from lxml import etree
from genshi.core import Markup
import readers

def target_decoy_generator(element_generator, decoy, ns):
    for ev,el in element_generator:
        if el.attrib['{%s}decoy' % ns['xmlns']] == decoy:
            strxml = etree.tostring(el)
            strxml = strxml.replace('xmlns="{0}" '.format(ns['xmlns']), '')
            strxml = strxml.replace('xmlns:p="{0}" '.format(ns['xmlns:p']), '')
            strxml = strxml.replace('xmlns:xsi="{0}" '.format(ns['xmlns:xsi']), '')
            strxml = Markup(strxml)
        else:
            continue
        el.clear()
        while el.getprevious() is not None:
            del el.getparent()[0]
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
    # first get all pep sequences
    pepseqs = {}
    pepgens = readers.get_peptides_multiple_fractions(input_files, ns)
    for pepgen in pepgens:
        for ac,el in pepgen:
            pepseqs[ el.attrib['{%s}peptide_id'] ] = 1

    sys.stdout.write('Filtering {0} unique peptide sequences')
    scores = {'q':'q_value', 'pep':'pep', 'p':'p_value', 'svm':'svm_score'}
    
    # now yield one unique peptide per seq. The looped over code could be in own
    # function to increase readability, but function calls are expensive.
    for pepseq in pepseqs:
        pepgens = readers.get_peptides_multiple_fractions(input_files, ns)
        # manually pop first peptide from generator to populate highest dict.
        ac,el = pepgen[0].next()
        featscore = float(el.xpath('xmlns:%s' % scores[score], namespaces=ns)[0].text)
        highest = {'score': featscore, 'pep_el': el}
        
        for pepgen in pepgens:
            for ac,el in pepgen:
                featscore = float(el.xpath('xmlns:%s' % scores[score], namespaces=ns)[0].text) 
                if score == 'svm': # greater than score is accepted
                    if featscore > highest['score']:
                        highest = {'pep_el': el, 'score': featscore}
                else: # lower than score is accepted
                    if featscore < highest['score']:
                        highest = {'pep_el': el, 'score': featscore}
        
        yield highest['pep_el']
