import sys

from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata


def merge_mzidtsvs(fns, header):
    for fn in fns:
        if header != tsvreader.get_tsv_header(fn):
            raise RuntimeError('Headers of TSV files to concatenate are '
                               'not identical')
    for psm in tsvreader.generate_tsv_lines_multifile(fns, header):
        yield psm


def get_splitfield(header, splitcol):
    """Returns column nr on which to split PSM table. Chooses from flags
    given via bioset and splitcol"""
    if splitcol == 'bioset':
        return mzidtsvdata.HEADER_SETNAME
    elif splitcol == 'TD':
        return mzidtsvdata.HEADER_TARGETDECOY
    else:
        try:
            splitcol = int(splitcol)
        except ValueError:
            print('ERROR: --splitcol must be an integer or "TD", or "bioset"')
            sys.exit(1)
        else:
            return header[splitcol - 1]


def generate_psms_split(psms, splitfield):
    """Loops PSMs and outputs dictionaries passed to writer. Dictionaries
    contain the PSMs and info to which split pool the
    respective PSM belongs"""
    for psm in psms:
        yield {'psm': psm,
               'split_pool': psm[splitfield],
               }
