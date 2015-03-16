from app.readers import tsv as readers
from app.dataformats import mzidtsv as mzidtsvdata


def create_header(oldheader, spectracol):
    newheadings = mzidtsvdata.SPECDATA_HEADER
    header = oldheader[:spectracol] + newheadings + oldheader[spectracol:]
    return header


def generate_psms_spectradata(lookup, tsvfn, oldheader,	spec_column):
    if spec_column is not None:
        speccol_head = oldheader[spec_column - 1]
    else:
        speccol_head = mzidtsvdata.HEADER_SPECFILE
    mzmlmap = lookup.get_mzmlfile_map()
    for psm in readers.generate_tsv_psms(tsvfn, oldheader):
        outpsm = {x: y for x, y in psm.items()}
        specfile = outpsm[speccol_head]
        scannr = outpsm[mzidtsvdata.HEADER_SCANNR]
        outpsm.update(lookup_spectra(lookup, mzmlmap[specfile], scannr))
        yield outpsm


def lookup_spectra(lookup, spectrafile_id, scannr):
    """Outputs dict with keys == spectradataheadernames,
    values == spectra data."""
    specdata = lookup.get_spectradata(spectrafile_id, scannr)
    return {mzidtsvdata.HEADER_SETNAME: specdata[0],
            mzidtsvdata.HEADER_RETENTION_TIME: specdata[1]}
