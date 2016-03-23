from app.actions.mzidtsv import percolator as prep
from app.drivers.mzidtsv import MzidTSVDriver
from app.drivers.options import mzidtsv_options


class MzidPercoTSVDriver(MzidTSVDriver):
    """
    Adds percolator data from mzid file to table.
    """
    outsuffix = '_percolated.tsv'
    command = 'percolator'
    commandhelp = ('Add percolator data to a  TSV with MSGF+ output. ')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['mzidfn'], mzidtsv_options))

    def get_psms(self):
        self.multipsm_per_scan = False
        if self.multipsm_per_scan:
            # FIXME not supported yet
            # Create mzid PSM/sequence sqlite (fn, scan, rank, sequence)
            pass
        self.header = prep.get_header_with_percolator(self.oldheader,
                                                      self.multipsm_per_scan)
        self.psms = prep.add_percolator_to_mzidtsv(self.mzidfn,
                                                   self.fn,
                                                   self.multipsm_per_scan,
                                                   self.oldheader)
