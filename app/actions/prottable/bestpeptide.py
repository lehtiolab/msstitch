from app.dataformats import peptable as peptabledata
from app.dataformats import prottable as prottabledata
from app.readers import tsv as reader
from app.actions.peptable.base import evaluate_peptide
from math import log


def generate_proteins(pepfn, proteins, pepheader, scorecol, minlog,
                      higherbetter=True):
    """Best peptide for each protein in a table"""
    protein_peptides = {}
    if minlog:
        higherbetter = False
    for psm in reader.generate_tsv_psms(pepfn, pepheader):
        p_acc = psm[peptabledata.HEADER_MASTERPROTEINS]
        if ';' in p_acc:
            continue
        protein_peptides = evaluate_peptide(protein_peptides, psm, p_acc,
                                            higherbetter, scorecol,
                                            fncol=False)
    if minlog:
        nextbestscore = min([pep['score'] for pep in protein_peptides.values()
                             if pep['score'] > 0])
        nextbestscore = -log(nextbestscore, 10)
    for protein in proteins:
        peptide = protein_peptides[protein[prottabledata.HEADER_PROTEIN]]
        if minlog:
            peptide['score'] = log_score(peptide['score'], nextbestscore)
        protein[prottabledata.HEADER_QSCORE] = str(
            peptide['score'])
        yield protein


def log_score(score, nextbestscore):
    try:
        logged_score = -log(score, 10)
    except ValueError:
        logged_score = nextbestscore + 2
    return logged_score
