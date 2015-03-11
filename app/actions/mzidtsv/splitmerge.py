import re

from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata


def merge_mzidtsvs(fns, header):
    for fn in fns:
        if header != tsvreader.get_tsv_header(fn):
            raise RuntimeError('Headers of TSV files to concatenate are '
                               'not identical')
    for psm in tsvreader.generate_tsv_lines_multifile(fns, header):
        yield psm


def get_splitcolnr(header, bioset, splitcol):
    """Returns column nr on which to split PSM table. Chooses from flags
    given via bioset and splitcol"""
    if bioset is not None:
        return header[mzidtsvdata.HEADER_SETNAME]
    elif splitcol is not None:
        return splitcol - 1
    else:
        raise RuntimeError('Must specify either --bioset or --splitcol')


def get_splitheader(oldheader, rename_cols=None, renamepattern=None):
    """Returns header for split files, adds a formatting string to
    certain columns if user has asked for this"""
    header = oldheader[:]
    if rename_cols is not None:
        for colnr in rename_cols:
            header[colnr - 1] = '{0}_' + header[colnr - 1]
    if renamepattern is not None:
        for colnr, col in enumerate(oldheader):
            if re.match('^{0}'.format(renamepattern), col):
                header[colnr] = '{0}_' + header[colnr]
    return header


def generate_psms_split(fn, oldheader, baseheader, bioset, splitcol):
    """Loops PSMs and outputs dictionaries passed to writer. Dictionaries
    contain the PSMs and info to which split pool the
    respective PSM belongs"""
    splitcolnr = get_splitcolnr(oldheader, bioset, splitcol)
    for psm in tsvreader.generate_tsv_psms(fn, baseheader):
        yield {'psm': psm,
               'split_pool':  ''.join(x for x in psm[oldheader[splitcolnr]]
                                      if x.isalnum() or x in '_.-')
               }
