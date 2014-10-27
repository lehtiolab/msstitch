import re
import itertools
from app.dataformats import mzidtsv as mzidtsvdata


def get_tsv_header(tsvfn):
    with open(tsvfn) as fp:
        return next(fp).strip().split('\t')


def generate_tsv_lines_multifile(fns, header):
    return itertools.chain.from_iterable([generate_tsv_psms(fn, header)
                                          for fn in fns])


def generate_tsv_psms(fn, header):
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
    """From a line, return list of proteins reported by Mzid2TSV. The line
    should not be unrolled."""
    proteins = line[mzidtsvdata.HEADER_PROTEIN].split(';')
    outproteins = []
    for protein in proteins:
        try:
            outproteins.append(protein[:protein.index('(')].strip())
        except ValueError:
            outproteins.append(protein)
    return outproteins


def get_pepproteins(line, unroll=False):
    """Returns from a PSM line peptide sequence,
    and other information about the PSM.
    Return values:
        specfn          -   str
        scan            -   str
        peptideseq      -   str
        score		-   str
        proteins        -   list of str
    """
    specfn = line[mzidtsvdata.HEADER_SPECFILE]
    scan = line[mzidtsvdata.HEADER_SCANNR]
    score = line[mzidtsvdata.HEADER_MSGFSCORE]
    peptideseq = line[mzidtsvdata.HEADER_PEPTIDE]
    if unroll and '.' in peptideseq:
        peptideseq = peptideseq.split('.')[1]
    proteins = get_proteins_from_psm(line)
    return specfn, scan, peptideseq, score, proteins


def strip_modifications(seq):
    nomodseq = re.sub('\d+', '', seq)
    nomodseq = re.sub('\+\.', '', seq)
    return nomodseq
