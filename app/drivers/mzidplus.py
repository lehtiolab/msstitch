from app.drivers import base
from app.preparation import mzidplus as prep
from app.writers import mzidplus as writers
from app.readers import basereader


class MzidPlusDriver(base.BaseDriver):
    def run(self):
        self.get_psms()
        self.write()

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_mzid_tsv(self.header, self.psms, outfn)


class MzidTSVConcatenateDriver(MzidPlusDriver):
    """Concatenates TSVs"""
    outsuffix = '_concat.tsv'

    def __init__(self, **kwargs):
        super(MzidTSVConcatenateDriver, self).__init__(**kwargs)
        self.allinfiles = [self.fn]
        self.allinfiles.extend(kwargs.get('multifile_input', None))

    def get_psms(self):
        self.header = prep.get_header_from_mzidtsv(self.fn)
        self.psms = basereader.generate_tsv_lines_multifile(self.allinfiles)


class MzidPercoTSVDriver(MzidPlusDriver):
    """
    Adds percolator data from mzid file to table.
    """
    outsuffix = '_percolated.tsv'

    def __init__(self, **kwargs):
        super(MzidPercoTSVDriver, self).__init__(**kwargs)
        self.idfn = kwargs.get('mzid', None)
        self.multipsm_per_scan = kwargs.get('allpsms', False)
        assert self.idfn is not None

    def get_psms(self):
        if self.multipsm_per_scan is True:
            # FIXME not supported yet
            # Create mzid PSM/sequence sqlite (fn, scan, rank, sequence)
            pass
        else:
            seqlookup = None

        self.header = prep.get_header_with_percolator(self.fn,
                                                      self.multipsm_per_scan)
        self.psms = prep.add_percolator_to_mzidtsv(self.idfn,
                                                   self.fn,
                                                   self.multipsm_per_scan,
                                                   self.header,
                                                   seqlookup)
