from app.dataformats import peptable as peptabledata
from app.dataformats import mzidtsv as psmtsvdata


def add_peptide(allpeps, psm, scorecol=False, fncol=None, new=False,
                track_psms=True):
    peptide = {'score': psm[scorecol],
               'line': psm,
               'psms': []
               }
    if track_psms:
        if not new:
            peptide['psms'] = allpeps[psm[peptabledata.HEADER_PEPTIDE]]['psms']
        peptide['psms'].append('{0}_{1}'.format(psm[fncol],
                                                psm[psmtsvdata.HEADER_SCANNR]))
    allpeps[psm[peptabledata.HEADER_PEPTIDE]] = peptide


def evaluate_peptide(peptides, psm, key, higherbetter, scorecol, fncol=None,
                     track_psms=True):
    try:
        existing_score = peptides[key]['score']
    except KeyError:
        add_peptide(peptides, psm, scorecol, fncol, True, track_psms)
    else:
        if higherbetter and psm[scorecol] > existing_score:
            add_peptide(peptides, psm, scorecol, fncol, track_psms=track_psms)
        elif not higherbetter and psm[scorecol] < existing_score:
            add_peptide(peptides, psm, scorecol, fncol, track_psms=track_psms)
    return peptides
