import re

from app.dataformats import mzidtsv as header
from app.actions.percolator.filters import sequence_to_filterseqs


def filter_psms_conf(psms, confkey, conflvl, lower_is_better):
    threshold = float(conflvl)
    for psm in psms:
        try:
            confval = float(psm[confkey])
        except (TypeError, ValueError):
            pass
        else:
            lower = confval < threshold
            if lower == lower_is_better:
                yield psm


def filter_psms_remove_set(psms, setnames):
    biosets = set(setnames)
    for psm in psms:
        if psm[header.HEADER_SETNAME] not in biosets:
            yield psm


def filter_known_searchspace(psms, lookup, seqcol, ntermwildcards, deamidation):
    """Yields peptides from generator as long as their sequence is not found in
    known search space dict. Useful for excluding peptides that are found in
    e.g. ENSEMBL or similar"""
    for psm in psms:
        seq_is_known = False
        seq = re.sub('[^A-Z]', '', psm[seqcol].upper())
        for seq in  sequence_to_filterseqs(seq, deamidation):
            if lookup.check_seq_exists(seq, ntermwildcards):
                seq_is_known = True
                break
        if not seq_is_known:
            yield psm


def filter_whole_proteins(psms, lookup, seqcol, deamidation, minpeplen, enforce_tryp):
    for psm in psms:
        seq_matches_protein = False
        seq = re.sub('[^A-Z]', '', psm[seqcol].upper())
        filterseqs = sequence_to_filterseqs(seq, deamidation)
        found_prots = lookup.get_proteins_from_peps(filterseqs, minpeplen)
        for pepseq, proteins in found_prots.items():
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
        if not seq_matches_protein:
            yield psm

