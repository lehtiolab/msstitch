import re
import os
import itertools
from app.dataformats import mzidtsv as mzidtsvdata
from app.dataformats import prottable as prottabledata


def get_tsv_header(tsvfn):
    with open(tsvfn) as fp:
        return next(fp).strip().split('\t')


def generate_tsv_lines_multifile(fns, header):
    return itertools.chain.from_iterable([generate_tsv_psms(fn, header)
                                          for fn in fns])


def generate_tsv_proteins(fn, header):
    return generate_split_tsv_lines(fn, header)


def generate_tsv_pep_protein_quants(fns):
    """Unlike generate_tsv_lines_multifile, this generates tsv lines
    from multiple files that may have different headers. Yields
    fn, header as well as quant data for each protein quant"""
    for fn in fns:
        header = get_tsv_header(fn)
        for pquant in generate_split_tsv_lines(fn, header):
            yield os.path.basename(fn), header, pquant


def generate_kronik_feats(fn):
    """Generates features from a Kronik output file"""
    header = get_tsv_header(fn)
    return generate_split_tsv_lines(fn, header)


def mzmlfn_kronikfeature_generator(mzmlfns, kronikfns):
    """Generates tuples of spectra filename and corresponding output
    features from kronik"""
    for mzmlfn, kronikfn in zip(mzmlfns, kronikfns):
        for quant_el in generate_kronik_feats(kronikfn):
            yield os.path.basename(mzmlfn), quant_el


def generate_tsv_psms(fn, header):
    return generate_split_tsv_lines(fn, header)


def generate_tsv_peptides(fn):
    header = get_tsv_header(fn)
    return generate_split_tsv_lines(fn, header)


def generate_split_tsv_lines(fn, header):
    """Returns dicts with header-keys and psm statistic values"""
    for line in generate_tsv_psms_line(fn):
        yield {x: y for (x, y) in zip(header, line.strip().split('\t'))}


def generate_tsv_psms_line(fn):
    with open(fn) as fp:
        next(fp)  # skip header
        for line in fp:
            yield line


def get_psm_id(line):
    return '{0}_{1}_{2}'.format(line[mzidtsvdata.HEADER_SPECFILE],
                                line[mzidtsvdata.HEADER_SCANNR],
                                line[mzidtsvdata.HEADER_PEPTIDE])


def get_proteins_from_psm(line):
    """From a line, return list of proteins reported by Mzid2TSV. When unrolled
    lines are given, this returns the single protein from the line."""
    proteins = line[mzidtsvdata.HEADER_PROTEIN].split(';')
    outproteins = []
    for protein in proteins:
        prepost_protein = re.sub('\(pre=.*post=.*\)', '', protein).strip()
        outproteins.append(prepost_protein)
    return outproteins


def get_psm(line, unroll=False, specfncol=None):
    """Returns from a PSM line peptide sequence,
    and other information about the PSM.
    Return values:
        specfn          -   str
        psm_id 	 	-   str
        scan            -   str
        peptideseq      -   str
        score		-   str
    """
    if specfncol is None:
        specfncol = mzidtsvdata.HEADER_SPECFILE
    specfn = line[specfncol]
    psm_id = get_psm_id(line)
    scan = line[mzidtsvdata.HEADER_SCANNR]
    peptideseq = get_psm_sequence(line, unroll)
    score = line[mzidtsvdata.HEADER_MSGFSCORE]
    return specfn, psm_id, scan, peptideseq, score


def get_psm_sequence(line, unroll=False):
    peptideseq = line[mzidtsvdata.HEADER_PEPTIDE]
    if unroll and '.' in peptideseq:
        peptideseq = peptideseq.split('.')[1]
    return peptideseq


def get_pepproteins(line):
    """Returns from a PSM line peptide sequence,
    and other information about the PSM.
    Return values:
        psm_id          -   str
        proteins        -   list of str
    """
    psm_id = get_psm_id(line)
    proteins = get_proteins_from_psm(line)
    return psm_id, proteins


def strip_modifications(seq):
    return re.sub('[-+.\d]', '', seq)


def get_cols_in_file(pattern, header, single_col=False):
    if pattern is None:
        return False
    cols_found = get_columns_by_pattern(header, pattern)
    if single_col:
        cols_found = cols_found[0]
    return cols_found


def get_columns_by_pattern(header, pattern):
    columns = []
    for field in header:
        if re.search(pattern, field) is not None:
            columns.append(field)
    if not columns:
        raise RuntimeError('Could not find fieldname in header with '
                           'pattern: {}'.format(pattern))
    return columns
