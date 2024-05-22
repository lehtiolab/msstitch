import re
import os
import itertools
from app.dataformats import prottable as prottabledata


def get_tsv_header(tsvfn):
    with open(tsvfn) as fp:
        return next(fp).strip('\n').split('\t')


def generate_tsv_lines_multifile(fns, header):
    return itertools.chain.from_iterable([generate_split_tsv_lines(fn, header)
                                          for fn in fns])


def generate_tsv_pep_protein_quants(fns):
    """Unlike generate_tsv_lines_multifile, this generates tsv lines
    from multiple files that may have different headers. Yields
    fn, header as well as quant data for each protein quant"""
    for fn in fns:
        basefn = os.path.basename(fn)
        header = get_tsv_header(fn)
        for pquant in generate_split_tsv_lines(fn, header):
            yield basefn, header, pquant


def generate_ms1_feats(fn):
    """Generates features from a Kronik output file"""
    header = get_tsv_header(fn)
    return generate_split_tsv_lines(fn, header)


def mzmlfn_tsvfeature_generator(mzmlfns, ms1fns):
    """Generates tuples of spectra filename and corresponding output
    features from dinosaur/kronik"""
    for mzmlfn, ms1fn in zip(mzmlfns, ms1fns):
        basefn = os.path.basename(mzmlfn)
        for quant_el in generate_ms1_feats(ms1fn):
            yield basefn, quant_el


def generate_split_tsv_lines(fn, header):
    """Returns dicts with header-keys and psm statistic values"""
    with open(fn) as fp:
        next(fp)  # skip header
        for line in fp:
            yield {x: y for (x, y) in zip(header, line.strip('\n').split('\t'))}


def strip_modifications(seq):
    # FIXME make method of PSM class??? maybe not
    return re.sub('[^A-Za-z]', '', seq)


def get_cols_in_file(pattern, header, single_col=False):
    if pattern is None:
        return False
    cols_found = get_columns_by_pattern(header, pattern)
    if single_col:
        cols_found = cols_found[0]
    return cols_found


def get_columns_by_combined_patterns(header, patterns):
    columns = []
    for field in header:
        if all(re.search(p, field) for p in patterns):
            columns.append(field)
    if not columns:
        raise RuntimeError('Could not find fieldname in header with '
                           f'combined patterns: {patterns}')
    return columns


def get_columns_by_pattern(header, pattern):
    columns = []
    for field in header:
        if re.search(pattern, field) is not None:
            columns.append(field)
    if not columns:
        raise RuntimeError('Could not find fieldname in header with '
                           'pattern: {}'.format(pattern))
    return columns
