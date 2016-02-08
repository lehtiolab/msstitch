from app.readers import mzidplus as readers
from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata


def add_percolator_to_mzidtsv(mzidfn, tsvfn, multipsm, oldheader):
    """Takes a MSGF+ tsv and corresponding mzId, adds percolatordata
    to tsv lines. Generator yields the lines. Multiple PSMs per scan
    can be delivered, in which case rank is also reported.
    """
    namespace = readers.get_mzid_namespace(mzidfn)
    try:
        xmlns = '{%s}' % namespace['xmlns']
    except TypeError:
        xmlns = ''
    specfnids = readers.get_mzid_specfile_ids(mzidfn, namespace)
    mzidpepmap = {}
    for peptide in readers.generate_mzid_peptides(mzidfn, namespace):
        pep_id, seq = readers.get_mzid_peptidedata(peptide, xmlns)
        mzidpepmap[pep_id] = seq
    mzidpercomap = {}
    for specid_data in readers.generate_mzid_spec_id_items(mzidfn, namespace,
                                                           xmlns, specfnids):
        scan, fn, pepid, spec_id = specid_data
        percodata = readers.get_specidentitem_percolator_data(spec_id, xmlns)
        try:
            mzidpercomap[fn][scan][mzidpepmap[pepid]] = percodata
        except KeyError:
            try:
                mzidpercomap[fn][scan] = {mzidpepmap[pepid]: percodata}
            except KeyError:
                mzidpercomap[fn] = {scan: {mzidpepmap[pepid]: percodata}}
    for line in tsvreader.generate_tsv_psms(tsvfn, oldheader):
        outline = {k: v for k, v in line.items()}
        fn = line[mzidtsvdata.HEADER_SPECFILE]
        scan = line[mzidtsvdata.HEADER_SCANNR]
        seq = line[mzidtsvdata.HEADER_PEPTIDE]
        outline.update(mzidpercomap[fn][scan][seq])
        yield outline


def get_header_with_percolator(oldheader, multipsm=False):
    percoheader = mzidtsvdata.PERCO_HEADER
    if multipsm is True:
        # FIXME should this be here???
        percoheader.append('rank')
    ix = oldheader.index(mzidtsvdata.HEADER_EVALUE) + 1
    return oldheader[:ix] + percoheader + oldheader[ix:]
