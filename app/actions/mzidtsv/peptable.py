from app.dataformats import mzidtsv as mzidtsvdata
from app.readers import tsv as reader


def get_peptable_header(oldheader):
    peptable_header = [mzidtsvdata.HEADER_LINKED_PSMS]
    ix = oldheader.index(mzidtsvdata.HEADER_PEPTIDE)
    return oldheader[:ix] + peptable_header + oldheader[ix:]


def generate_peptides(tsvfn, oldheader, scorecol, fncol=False,
                      higherbetter=True):
    if fncol is None:
        fncol = mzidtsvdata.HEADER_SPECFILE
    peptides = {}
    for psm in reader.generate_tsv_psms(tsvfn, oldheader):
        try:
            existing_score = peptides[psm[mzidtsvdata.HEADER_PEPTIDE]]['score']
        except KeyError:
            add_peptide(peptides, psm, fncol, scorecol,
                        new=True)
        else:
            if higherbetter and existing_score < psm[scorecol]:
                add_peptide(peptides, psm, fncol, scorecol)
    for peptide in peptides.values():
        peptide['line'][mzidtsvdata.HEADER_LINKED_PSMS] = '; '.join(
            peptide['psms'])
        yield peptide['line']


def add_peptide(allpeps, psm, fncol, scorecol=False, new=False):
    peptide = {'score': psm[scorecol],
               'line': psm,
               'psms': []
               }
    if not new:
        peptide['psms'] = allpeps[psm[mzidtsvdata.HEADER_PEPTIDE]]['psms']
    peptide['psms'].append('{0}_{1}'.format(psm[fncol],
                                            psm[mzidtsvdata.HEADER_SCANNR]))
    allpeps[psm[mzidtsvdata.HEADER_PEPTIDE]] = peptide
    return allpeps
