import itertools
from hashlib import md5
TSV_SPECFN_COL = 1
TSV_SCAN_COL = 2
TSV_PEPTIDE_COL = 8
TSV_PROTEIN_COL = 9


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


def get_mzidtsv_lines_scannr_specfn(fn):
    """Returns generator of lines of tsv, skipping header, as split lists,
    and a tuple containing (spectra file, scan nr)."""
    header = get_tsv_header(fn)
    assert header[TSV_SCAN_COL] in ['ScanNum', 'scannr', 'scan_nr',
                                    'Scan number']
    assert header[TSV_SPECFN_COL] in ['#SpecFile', 'spectra_file', 'specfile']
    for line in generate_tsv_psms_line(fn):
        line = line.strip().split('\t')
        yield line, (line[TSV_SPECFN_COL], line[TSV_SCAN_COL])


def get_multiple_proteins(line):
    """From a line, return list of proteins reported by Mzid2TSV. The line
    should not be unrolled."""
    proteins = line[TSV_PROTEIN_COL].split(';')
    return [x[:x.index('(')].strip() for x in proteins]


def get_unrolled_proteins(line):
    """From a line, return the protein reported by Mzid2TSV as a list. This
    function applies to tsvs that have one protein per line (unrolled)"""
    return [line[TSV_PROTEIN_COL]]




def get_peptide_proteins(line, specfn, scannr, unroll=False):
    """From a line, generate a peptide_id (MD5 of specfile and scannr).
    Returns that id with the peptide sequence and protein accessions.
    Return values:
        peptide_id      -   str
        peptideseq      -   str
        proteins        -   list of str
    """
    peptideseq = line[TSV_PEPTIDE_COL]
    peptide_id = md5('{0}{1}'.format(specfn, scannr)
                     .encode('utf-8')).hexdigest()
    if unroll:
        peptideseq = peptideseq.split('.')[1]
        proteins = get_unrolled_proteins(line)
    else:
        proteins = get_multiple_proteins(line)
    return peptide_id, peptideseq, proteins
