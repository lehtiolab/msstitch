import app.preparation.mzidtsv.percolator as prep
from app.drivers.mzidtsv import MzidTSVDriver


class MzidPercoTSVDriver(MzidTSVDriver):
    """
    Adds percolator data from mzid file to table.
    """
    outsuffix = '_percolated.tsv'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
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

        self.header = prep.get_header_with_percolator(self.oldheader,
                                                      self.multipsm_per_scan)
        self.psms = prep.add_percolator_to_mzidtsv(self.idfn,
                                                   self.fn,
                                                   self.multipsm_per_scan,
                                                   self.oldheader,
                                                   self.header,
                                                   seqlookup)
