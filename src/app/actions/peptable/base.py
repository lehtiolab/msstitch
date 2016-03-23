from app.dataformats import mzidtsv as psmtsvdata


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
                                            psm[psmtsvdata.HEADER_SCANNR]))


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
