from numpy import polyfit
from math import log

from app.readers import tsv as reader
from app.dataformats import peptable as peptabledata


def recalculate_qvals_linear_model(fn, scorecol, qvalcol, qvalthreshold):
    slope, intercept = fit_linear_model(fn, scorecol, qvalcol, qvalthreshold)
    for peptide in reader.generate_tsv_peptides(fn):
        outpeptide = {k: v for k, v in peptide.items()}
        score = float(outpeptide[scorecol])
        if not slope:
            qval = 'NA'
        else:
            qval = 10 ** (slope * score + intercept)
        outpeptide[peptabledata.HEADER_QVAL_MODELED] = str(qval)
        yield outpeptide


def fit_linear_model(fn, scorecol, qvalcol, qvalthreshold):
    pepq = []
    pepscore = []
    for peptide in reader.generate_tsv_peptides(fn):
        qval = float(peptide[qvalcol])
        if qval > qvalthreshold:
            pepq.append(log(qval, 10))
            pepscore.append(float(peptide[scorecol]))
    if len(pepq) < 10:
        slope, intercept = False, False
        print('Could not fit linear model through q-values, as only {} q-values '
                'was/were above the q-value threshold of {} for inclusion'.format(len(pepq), qvalthreshold))
    else:
        slope, intercept = polyfit(pepscore, pepq, deg=1)
        print('Fitted linear model through qvalues (above {}) vs score. '
              'Slope={}, intercept={}'.format(qvalthreshold, slope, intercept))
    return slope, intercept
