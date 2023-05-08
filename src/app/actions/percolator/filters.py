import re
from itertools import product

from app.readers import percolator as reader
from app.readers import xmlformatting as formatting


def filter_whole_proteins(elements, lookup, seqtype, ns, deamidation, minpeplen, enforce_tryp):
    # Remove duplicate sequences by using them as keys
    for element in elements:
        seq_matches_protein = False
        element_seqs = get_seqs_from_element(element, seqtype, ns, deamidation)
        element_prots = lookup.get_proteins_from_peps(element_seqs, minpeplen)
        for pepseq, proteins in element_prots.items():
            for prot_id, pos, protseq in proteins:
                if pepseq in protseq.replace('L', 'I'):
                    if enforce_tryp and (pos == 0 or not set(
                            [pepseq[-1],
                             protseq[pos - 1]]).difference(['K', 'R'])):
                        # pepseq is tryptic on both ends, or
                        # pepseq is an N-term peptide),
                        # matches to protein seq so remove
                        seq_matches_protein = True
                        break
                    elif not enforce_tryp:
                        seq_matches_protein = True
                        break
        if seq_matches_protein:
            formatting.clear_el(element)
        else:
            yield formatting.string_and_clear(element, ns)


def filter_known_searchspace(elements, seqtype, lookup, ns, ntermwildcards,
                             deamidation):
    """Yields peptides from generator as long as their sequence is not found in
    known search space dict. Useful for excluding peptides that are found in
    e.g. ENSEMBL or similar"""
    for element in elements:
        seq_is_known = False
        for seq in get_seqs_from_element(element, seqtype, ns, deamidation):
            if lookup.check_seq_exists(seq, ntermwildcards):
                seq_is_known = True
                break
        if seq_is_known:
            formatting.clear_el(element)
        else:
            yield formatting.string_and_clear(element, ns)


def sequence_to_filterseqs(seq, deamidation):
    # Exchange leucines for isoleucines since MS can't differ and we
    # don't want to find 'novel' peptides which only have a difference
    # in this amino acid
    seq = seq.replace('L', 'I')
    if deamidation:
        options = [(c,) if c != 'D' else ('D', 'N') for c in seq]
        return list(''.join(o) for o in product(*options))
    else:
        return [seq]


def get_seqs_from_element(element, seqtype, ns, deamidation):
    seq = {'psm': reader.get_psm_seq, 'pep': reader.get_peptide_seq}[seqtype](element, ns).upper()
    seq = re.sub('\[UNIMOD:\d*\]', '', seq)
    return sequence_to_filterseqs(seq, deamidation)
