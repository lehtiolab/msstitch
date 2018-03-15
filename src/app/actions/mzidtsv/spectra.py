from app.readers import tsv as readers
from app.dataformats import mzidtsv as mzidtsvdata


def create_header(oldheader, spectracol):
    newheadings = mzidtsvdata.SPECDATA_HEADER
    header = oldheader[:spectracol] + newheadings + oldheader[spectracol:]
    return header


def generate_psms_spectradata(lookup, tsvfn, oldheader):
    psm_specdata = zip(enumerate(readers.generate_tsv_psms(tsvfn, oldheader)),
                       lookup.get_exp_spectra_data_rows())
    for (row, psm), specdata in psm_specdata:
        outpsm = {x: y for x, y in psm.items()}
        if row == int(specdata[0]):
            outpsm.update({mzidtsvdata.HEADER_SETNAME: specdata[1],
                           mzidtsvdata.HEADER_RETENTION_TIME: str(specdata[2]),
                           mzidtsvdata.HEADER_INJECTION_TIME: str(specdata[3]),
                           })
        else:
            raise RuntimeError('PSM with row nr {} has no rownr in DB. '
                               'Current DB row is {}'.format(row, specdata[0]))
        yield outpsm
