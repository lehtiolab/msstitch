from app.actions.mzidtsv import filter_confidence as prep
from app.drivers.mzidtsv import MzidTSVDriver


class ConfidenceFilterDriver(MzidTSVDriver):
    outsuffix = '_filtconf.txt'
    command = 'conffilt'

    def get_psms(self):
        confkey = self.oldheader[int(self.confcol) - 1]
        self.header = self.oldheader[:]
        self.psms = prep.generate_psms(self.fn,
                                       self.oldheader,
                                       confkey,
                                       self.conflvl,
                                       self.lowerbetter)
