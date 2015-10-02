from app.dataformats import mzidtsv as psmtsvdata


def add_peptide(allpeps, psm, key, score, fncol=None, new=False,
                track_psms=True):
    peptide = {'score': score,
               'line': psm,
               'psms': []
               }
    if track_psms:
        if not new:
            peptide['psms'] = allpeps[key]['psms']
        peptide['psms'].append('{0}_{1}'.format(psm[fncol],
                                                psm[psmtsvdata.HEADER_SCANNR]))
    allpeps[key] = peptide


def evaluate_peptide(peptides, psm, key, higherbetter, scorecol, fncol=None,
                     track_psms=True):
    try:
        score = float(psm[scorecol])
    except ValueError:
        # If score is NA or similar, dont use this PSM 
        return peptides
    try:
        existing_score = peptides[key]['score']
    except KeyError:
        add_peptide(peptides, psm, key, score, fncol, True, track_psms)
    else:
        if higherbetter and score > existing_score:
            add_peptide(peptides, psm, key, score, fncol,
                        track_psms=track_psms)
        elif not higherbetter and score < existing_score:
            add_peptide(peptides, psm, key, score, fncol,
                        track_psms=track_psms)
    return peptides
