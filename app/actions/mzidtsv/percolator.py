from app.readers import mzidplus as readers
from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata


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
                              oldheader, seqdb=None):
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
    #  not match. Also error at not found lines.
    specresult, specdata = get_specresult_data(specresults, specfnids)
    writelines = []
    for line in tsvreader.generate_tsv_psms(tsvfn, oldheader):
        while True:
            if line[mzidtsvdata.HEADER_SCANNR] == specdata['scan'] and \
               line[mzidtsvdata.HEADER_SPECFILE] == specdata['fn']:
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
        percoheader.append('rank')
    ix = oldheader.index(mzidtsvdata.HEADER_EVALUE) + 1
    return oldheader[:ix] + percoheader + oldheader[ix:]
