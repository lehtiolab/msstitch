from app.actions.mzidtsv import proteingrouping as prep
from app.drivers.mzidtsv import MzidTSVDriver


class ProteinGroupDriver(MzidTSVDriver):
    outsuffix = '_protgroups.txt'
    lookuptype = 'proteingroups'

    def get_psms(self):
        confkey = self.oldheader[int(self.confcol) - 1]
        self.header = prep.get_header_with_proteingroups(self.oldheader)
        self.psms = prep.generate_psms_with_proteingroups(self.fn,
                                                          self.oldheader,
                                                          self.header,
                                                          self.lookup,
                                                          confkey,
                                                          self.conflvl,
                                                          self.lowerbetter,
                                                          self.unroll)
