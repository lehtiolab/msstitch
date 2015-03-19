import re
from lxml import etree
from app import formatting
from app.lookups.sqlite import searchspace as sqlite
from app.readers import pycolator as reader


def get_either_seq(seqtype, element, ns):
    get_seq_map = {'psm': get_psm_seq,
                   'pep': get_peptide_seq,
                   }
    return get_seq_map[seqtype](element, ns)


def get_peptide_seq(peptide, ns):
    return peptide.attrib['{%s}peptide_id' % ns['xmlns']]


def get_psm_seq(psm, ns):
    return psm.find('{%s}peptide_seq' % ns['xmlns']).attrib['seq']


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


def filter_peptide_length(features, elementtype, ns, minlen=0, maxlen=None):
    minlen = int(minlen)
    if maxlen is None:
        maxlen = float('inf')
    else:
        maxlen = int(maxlen)
    for feat in features:
        seq = get_either_seq(elementtype, feat, ns)
        seq = strip_modifications(seq)
        if len(seq) > minlen and len(seq) < maxlen:
            yield formatting.string_and_clear(feat, ns)
        else:
            formatting.clear_el(feat)


def strip_modifications(seq):
    return re.sub('\[UNIMOD:\d*\]', '', seq)


def filter_known_searchspace(elements, seqtype, searchspace,
                             ns, ntermwildcards):
    """Yields peptides from generator as long as their sequence is not found in
    known search space dict. Useful for excluding peptides that are found in
    e.g. ENSEMBL or similar"""
    lookup = sqlite.SearchSpaceDB(searchspace)

    for element in elements:
        seq = get_either_seq(seqtype, element, ns)
        seq = strip_modifications(seq)
        # Exchange leucines for isoleucines since MS can't differ and we
        # don't want to find 'novel' peptides which only have a difference
        # in this amino acid
        seq = seq.replace('L', 'I')
        if not lookup.check_seq_exists(seq, ntermwildcards):
            yield formatting.string_and_clear(element, ns)
        else:
            formatting.clear_el(element)
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


def merge_peptides(fns, ns):
    """Loops peptides from multiple files, fetches PSMs from
    sequence:PSM map, outputs correctly PSM mapped peptides"""
    peptides_to_map = reader.generate_peptides_multiple_fractions(fns, ns)
    psmmap = create_merge_psm_map(peptides_to_map, ns)
    peptides = reader.generate_peptides_multiple_fractions(fns, ns)
    for peptide in peptides:
        seq = get_peptide_seq(peptide, ns)
        psm_ids = get_psm_ids_from_peptide(peptide, ns)
        # remove current psm ids, repopulate with stored ones
        psm_ids.clear()
        for new_psm_id in psmmap[seq]:
            etree.SubElement(psm_ids, 'psm_id').text = new_psm_id
        yield formatting.string_and_clear(peptide, ns)


def create_merge_psm_map(peptides, ns):
    """Loops through peptides, stores sequences mapped to PSM ids."""
    psmmap = {}
    for peptide in peptides:
        seq = get_peptide_seq(peptide, ns)
        psm_ids = get_psm_ids_from_peptide(peptide, ns)
        for psm_id in psm_ids:
            try:
                psmmap[seq][psm_id.text] = 1
            except KeyError:
                psmmap[seq] = {psm_id.text: 2}
    for seq, psm_id_dict in psmmap.items():
        psmmap[seq] = [x for x in psm_id_dict]
    return psmmap


def get_psm_ids_from_peptide(peptide, ns):
    return peptide.xpath('xmlns:psm_ids', namespaces=ns)[0]