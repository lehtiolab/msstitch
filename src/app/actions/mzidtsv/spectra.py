import re

from app.readers import tsv as readers
from app.dataformats import mzidtsv as mzidtsvdata


def create_header(oldheader, spectracol, bioset, miscleav):
    newheadings = mzidtsvdata.MOREDATA_HEADER
    if bioset:
        newheadings = [mzidtsvdata.HEADER_SETNAME] + newheadings
    if miscleav:
        newheadings.append(mzidtsvdata.HEADER_MISSED_CLEAVAGE)
    header = oldheader[:spectracol] + newheadings + oldheader[spectracol:]
    return header


def count_missed_cleavage(full_pepseq, count=0):
    '''Regex .*[KR][^P] matches until the end and checks if there is a final
    charachter so this will not match the tryptic residue'''
    pepseq = re.sub('[\+\-]\d*.\d*', '', full_pepseq)
    match = re.match('.*[KR][^P]', pepseq)
    if match:
        count += 1
        return count_missed_cleavage(match.group()[:-1], count)
    else:
        return count


def generate_psms_spectradata(lookup, tsvfn, oldheader, bioset, miscleav):
    psm_specdata = zip(enumerate(readers.generate_tsv_psms(tsvfn, oldheader)),
                       lookup.get_exp_spectra_data_rows())
    for (row, psm), specdata in psm_specdata:
        outpsm = {x: y for x, y in psm.items()}
        if row == int(specdata[0]):
            outpsm.update({mzidtsvdata.HEADER_RETENTION_TIME: str(specdata[2]),
                           mzidtsvdata.HEADER_INJECTION_TIME: str(specdata[3]),
                           })
            if bioset:
                outpsm[mzidtsvdata.HEADER_SETNAME] = str(specdata[1])
        else:
            raise RuntimeError('PSM with row nr {} has no rownr in DB. '
                               'Current DB row is {}'.format(row, specdata[0]))
        if miscleav:
            outpsm[mzidtsvdata.HEADER_MISSED_CLEAVAGE] = count_missed_cleavage(outpsm[mzidtsvdata.HEADER_PEPTIDE])
        yield outpsm
