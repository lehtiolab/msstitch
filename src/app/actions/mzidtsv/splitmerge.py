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
    if bioset:
        return header.index(mzidtsvdata.HEADER_SETNAME)
    elif splitcol is not None:
        return splitcol - 1
    else:
        raise RuntimeError('Must specify either --bioset or --splitcol')


def generate_psms_split(fn, oldheader, bioset, splitcol):
    """Loops PSMs and outputs dictionaries passed to writer. Dictionaries
    contain the PSMs and info to which split pool the
    respective PSM belongs"""
    try:
        splitcolnr = get_splitcolnr(oldheader, bioset, splitcol)
    except IndexError:
        raise RuntimeError('Cannot find bioset header column in '
                           'input file {}, though --bioset has '
                           'been passed'.format(fn))
    for psm in tsvreader.generate_tsv_psms(fn, oldheader):
        yield {'psm': psm,
               'split_pool': psm[oldheader[splitcolnr]]
               }
