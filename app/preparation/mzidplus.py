from app.readers import mzidplus as readers
from app.readers import tsv as tsvreader
HEADER_SPECFILE = '#SpecFile'
HEADER_SCANNR = 'ScanNum'
HEADER_PROTEIN = 'Protein'  # protein column as specified by Mzid2TSV
HEADER_MASTER_PROT = 'Master protein(s)'
HEADER_PG_CONTENT = 'Protein group(s) content'
HEADER_PG_AMOUNT_PROTEIN_HITS = 'Amount of matching proteins in group(s)'
HEADER_PG = [HEADER_MASTER_PROT, HEADER_PG_CONTENT,
             HEADER_PG_AMOUNT_PROTEIN_HITS]
HEADER_PEPQVAL = 'PepQValue'


def merge_mzidtsvs(fns, header):
    for fn in fns:
        if header != tsvreader.get_tsv_header(fn):
            raise RuntimeError('Headers of TSV files to concatenate are '
                               'not identical')
    for psm in tsvreader.generate_tsv_lines_multifile(fns, header):
        yield psm


def get_percoline(specresult, namespace, line, multipsm, seqdb):
    """Extracts percolator data from specresult and returns a dict"""
    out = line
    out.update({'rank': None})
    try:
        xmlns = '{%s}' % namespace['xmlns']
    except TypeError:
        xmlns = ''
    if multipsm is True:
        pass  # FIXME support later
        # loop through psms in specresult
        # check line sequence (without mods) in seqdb with psm
        # get percodata,
        # percoline = [line-with-correct-rank]
    else:  # only the first element
        perco = readers.get_specidentitem_percolator_data(
            specresult.find('{0}SpectrumIdentificationItem'.format(xmlns)),
            namespace)
    out.update(perco)
    return out


def get_specresult_data(specresults, id_fnlookup):
    specresult = next(specresults)
    scannr = readers.get_specresult_scan_nr(specresult)
    mzmlid = readers.get_specresult_mzml_id(specresult)
    return specresult, {'scan': scannr, 'fn': id_fnlookup[mzmlid]}


def add_percolator_to_mzidtsv(mzidfn, tsvfn, multipsm,
                              oldheader, header, seqdb=None):
    """Takes a MSGF+ tsv and corresponding mzId, adds percolatordata
    to tsv lines. Generator yields the lines. Multiple PSMs per scan
    can be delivered, in which case rank is also reported.
    """
    namespace = readers.get_mzid_namespace(mzidfn)
    specfnids = readers.get_mzid_specfile_ids(mzidfn, namespace)
    specresults = readers.mzid_spec_result_generator(mzidfn, namespace)
    # multiple lines can belong to one specresult, so we use a nested
    # for/while-true-break construction.
    # FIXME we assume best ranking is first line. Fix this in
    # FIXME get header names instead of positions!
    # FIXME we should count amounts of specresult/line and throw error if
    # they do not match. Also error at not found lines.
    specresult, specdata = get_specresult_data(specresults, specfnids)
    writelines = []
    for line in tsvreader.generate_tsv_psms(tsvfn, oldheader):
        while True:
            if line[HEADER_SCANNR] == specdata['scan'] and \
               line[HEADER_SPECFILE] == specdata['fn']:
                outline = get_percoline(specresult, namespace, line,
                                        multipsm, seqdb)
                writelines.append(outline)
                break
            else:
                specresult, specdata = get_specresult_data(specresults,
                                                           specfnids)
        for outline in writelines:
            yield outline
        writelines = []
    # write last lines
    for outline in writelines:
        yield outline


def get_header_with_percolator(oldheader, multipsm=False):
    percoheader = readers.PERCO_HEADER
    if multipsm is True:
        # FIXME should this be here???
        header.append('rank')
    header.extend(readers.PERCO_HEADER)
    return header
