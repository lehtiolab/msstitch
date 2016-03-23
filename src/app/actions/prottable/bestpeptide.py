from app.dataformats import peptable as peptabledata
from app.dataformats import prottable as prottabledata
from app.readers import tsv as reader
from app.actions.peptable.base import evaluate_peptide
from math import log


def generate_proteins(pepfn, proteins, pepheader, scorecol, minlog,
                      higherbetter=True, protcol=False):
    """Best peptide for each protein in a table"""
    protein_peptides = {}
    if minlog:
        higherbetter = False
    if not protcol:
        protcol = peptabledata.HEADER_MASTERPROTEINS
    for psm in reader.generate_tsv_psms(pepfn, pepheader):
        p_acc = psm[protcol]
        if ';' in p_acc:
            continue
        protein_peptides = evaluate_peptide(protein_peptides, psm, p_acc,
                                            higherbetter, scorecol,
                                            fncol=False)
    if minlog:
        try:
            nextbestscore = min([pep['score'] for pep in
                                 protein_peptides.values()
                                 if pep['score'] > 0])
        except ValueError:
            import sys
            sys.stderr.write('Cannot find score of type {} which is above 0. '
                             'Only scores above zero can have a -log value. '
                             'Exiting.'.format(scorecol))
            sys.exit(1)
        nextbestscore = -log(nextbestscore, 10)
    for protein in proteins:
        try:
            peptide = protein_peptides[protein[prottabledata.HEADER_PROTEIN]]
        except KeyError:
            print('WARNING - protein {} not found in peptide '
                  'table'.format(protein[prottabledata.HEADER_PROTEIN]))
            peptide = {'score': 'NA'}
        if minlog and peptide['score'] != 'NA':
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
