from statistics import median
from numpy import polyfit
from math import log
import sys

from app.dataformats import mzidtsv as mzidtsvdata
from app.dataformats import peptable as peptabledata
from app.readers import tsv as reader


def get_quantcols(pattern, oldheader, coltype):
    """Searches for quantification columns using pattern and header list.
    Calls reader function to do regexp. Returns a single column for
    precursor quant."""
    if pattern is None:
       return False
    if coltype == 'precur':
        return reader.get_cols_in_file(pattern, oldheader, single_col=True)


def generate_peptides(tsvfn, oldheader, switch_map, scorecol, precurquantcol,
        fncol=None, higherbetter=True):
    if fncol is None:
        fncol = mzidtsvdata.HEADER_SPECFILE
    peptides = {}
    for psm in reader.generate_split_tsv_lines(tsvfn, oldheader):
        for oldkey, newkey in switch_map.items():
            try:
                psm[newkey] = psm.pop(oldkey)
            except KeyError:
                pass
        pepseq = psm[peptabledata.HEADER_PEPTIDE]
        peptides = evaluate_peptide(peptides, psm, pepseq, higherbetter,
                                    scorecol, fncol)
        add_quant_values(peptides, psm, precurquantcol)
    for peptide in peptides.values():
        peptide['line'][peptabledata.HEADER_LINKED_PSMS] = '; '.join(
            peptide['psms'])
        for qtype, pepquant in peptide['quant'].items():
            peptide['line'].update(parse_quant_data(qtype, pepquant))
        yield peptide['line']


def parse_quant_data(qtype, pepquant):
    if qtype == 'precur':
        quants = {peptabledata.HEADER_AREA: get_peptide_quant(pepquant, 'precur')}
    return quants


def get_peptide_quant(quantdata, quanttype):
    """Parses lists of quantdata and returns maxvalue from them. Strips NA"""
    parsefnx = {'precur': max}
    quantfloats = []
    for q in quantdata:
        try:
            quantfloats.append(float(q))
        except(TypeError, ValueError):
            pass
    if not quantfloats:
        return 'NA'
    return str(parsefnx[quanttype](quantfloats))


def add_quant_values(allpeps, psm, precurq_field):
    try:
        quants = allpeps[psm[peptabledata.HEADER_PEPTIDE]]['quant']
    except KeyError:
        quants = {}
    if precurq_field:
        try:
            quants['precur'].append(psm[precurq_field])
        except KeyError:
            quants['precur'] = [psm[precurq_field]]
    try:
        allpeps[psm[peptabledata.HEADER_PEPTIDE]]['quant'] = quants
    except KeyError:
        # PSM sequence does not exist in allpeps because has no good-scoring peptides
        pass


def recalculate_qvals_linear_model(peptides, scorecol, qvalthreshold, minpeptidenr):
    peptides = [x for x in peptides]
    slope, intercept = fit_linear_model(peptides, scorecol, qvalthreshold,
            minpeptidenr)
    for peptide in peptides:
        outpeptide = {k: v for k, v in peptide.items()}
        score = float(outpeptide[scorecol])
        if not slope:
            qval = 'NA'
        else:
            qval = 10 ** (slope * score + intercept)
        outpeptide[peptabledata.HEADER_QVAL_MODELED] = str(qval)
        yield outpeptide


def fit_linear_model(peptides, scorecol, qvalthreshold, minpeptidenr):
    pepq = []
    pepscore = []
    for peptide in peptides:
        try:
            qval = float(peptide[peptabledata.HEADER_QVAL])
        except KeyError:
            print('PSM table has no {} column, cannot perform linear modeling '
                    'of q-values. Exiting.'.format(peptabledata.HEADER_QVAL))
            sys.exit(1)
        except ValueError:
            print('Q-value found which is not a number: {}, please check the '
                    'PSM table. Exiting.'.format(peptide[peptabledata.HEADER_QVAL]))
            sys.exit(1)
        if qval > qvalthreshold:
            pepq.append(log(qval, 10))
            pepscore.append(float(peptide[scorecol]))
    if len(pepq) < minpeptidenr:
        slope, intercept = False, False
        print('Could not fit linear model through q-values, as only {} q-values '
                '({} are needed) was/were above the q-value threshold of {} '
                'for modeling inclusion (model is created from highest '
                'q-values to approximate the stepped lower ones)'.format(
                    len(pepq), minpeptidenr, qvalthreshold))
    else:
        slope, intercept = polyfit(pepscore, pepq, deg=1)
        print('Fitted linear model through qvalues (above {}) vs score. '
              'Slope={}, intercept={}'.format(qvalthreshold, slope, intercept))
    return slope, intercept


def add_peptide(allpeps, psm, key, score, fncol=False, new=False):
    peptide = {'score': score,
               'line': psm,
               'psms': []
               }
    if fncol:
        if not new:
            peptide['psms'] = allpeps[key]['psms']
    allpeps[key] = peptide


def add_psm_link(peptide, psm, fncol):
    peptide['psms'].append('{0}_{1}'.format(psm[fncol],
                                            psm[mzidtsvdata.HEADER_SPECSCANID]))


def evaluate_peptide(peptides, psm, key, higherbetter, scorecol, fncol=False):
    try:
        score = float(psm[scorecol])
    except ValueError:
        # If score is NA or similar, dont use this PSM
        return peptides
    try:
        existing_score = peptides[key]['score']
    except KeyError:
        add_peptide(peptides, psm, key, score, fncol, new=True)
    else:
        if higherbetter and score > existing_score:
            add_peptide(peptides, psm, key, score, fncol)
        elif not higherbetter and score < existing_score:
            add_peptide(peptides, psm, key, score, fncol)
    if fncol:
        add_psm_link(peptides[key], psm, fncol)
    return peptides
