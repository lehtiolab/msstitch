import re
from app import formatting
from app import sqlite


def get_peptide_seq(peptide, ns):
    return peptide.attrib['{%s}peptide_id' % ns['xmlns']]


def target_decoy_generator(element_generator, decoy, ns):
    for el in element_generator:
        if el.attrib['{%s}decoy' % ns['xmlns']] == decoy:
            yield formatting.string_and_clear(el, ns)
        else:
            formatting.clear_el(el)


def split_target_decoy(elements, ns, filter_type):
    td = {'target': 'false', 'decoy': 'true'}
    feats_to_process = {'psm': None, 'peptide': None}
    for feat in feats_to_process:
        feats_to_process[feat] = target_decoy_generator(
            elements[feat], td[filter_type], ns)
    return feats_to_process


def get_score(elements, ns, scoretype='svm_score'):
    for el in elements:
        score = el.xpath('xmlns:{0}'.format(scoretype), namespaces=ns)[0].text
        formatting.clear_el(el)
        yield score


def filter_peptide_length(peptides, ns, minlen=0, maxlen=None):
    minlen = int(minlen)
    if maxlen is None:
        maxlen = float('inf')
    else:
        maxlen = int(maxlen)
    for pep in peptides:
        seq = pep.attrib['{%s}peptide_id' % ns['xmlns']]
        if len(seq) > minlen and len(seq) < maxlen:
            yield pep
        else:
            formatting.clear_el(pep)


def filter_known_searchspace(peptides, searchspace, ns, ntermwildcards):
    """Yields peptides from generator as long as their sequence is not found in
    known search space dict. Useful for excluding peptides that are found in
    e.g. ENSEMBL or similar"""
    lookup = sqlite.SearchSpaceDB(searchspace)

    for peptide in peptides:
        seq = get_peptide_seq(peptide, ns)
        # Exchange leucines for isoleucines since MS can't differ and we
        # don't want to find 'novel' peptides which only have a difference
        # in this amino acid
        seq = seq.replace('L', 'I')
        # Loose modifications
        seq = re.sub('\[UNIMOD:\d*\]', '', seq)
        if not lookup.check_seq_exists(seq, ntermwildcards):
            yield peptide
        else:
            formatting.clear_el(peptide)
    lookup.close_connection()


def filter_unique_peptides(peptides, score, ns):
    """ Filters unique peptides from multiple Percolator output XML files.
        Takes a dir with a set of XMLs, a score to filter on and a namespace.
        Outputs an ElementTree.
    """
    scores = {'q': 'q_value',
              'pep': 'pep',
              'p': 'p_value',
              'svm': 'svm_score'}
    highest = {}
    for el in peptides:
        featscore = float(el.xpath('xmlns:%s' % scores[score],
                                   namespaces=ns)[0].text)
        seq = get_peptide_seq(el, ns)

        if seq not in highest:
            highest[seq] = {
                'pep_el': formatting.stringify_strip_namespace_declaration(
                    el, ns), 'score': featscore}
        if score == 'svm':  # greater than score is accepted
            if featscore > highest[seq]['score']:
                highest[seq] = {
                    'pep_el':
                    formatting.stringify_strip_namespace_declaration(el, ns),
                    'score': featscore}
        else:  # lower than score is accepted
            if featscore < highest[seq]['score']:
                highest[seq] = {
                    'pep_el':
                    formatting.stringify_strip_namespace_declaration(el, ns),
                    'score': featscore}
        formatting.clear_el(el)

    for pep in list(highest.values()):
        yield pep['pep_el']
