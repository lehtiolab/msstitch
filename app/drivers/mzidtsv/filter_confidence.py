from app.actions.mzidtsv import filter_confidence as prep
from app.drivers.mzidtsv import MzidTSVDriver
from app.drivers.options import mzidtsv_options


class ConfidenceFilterDriver(MzidTSVDriver):
    outsuffix = '_filtconf.txt'
    command = 'conffilt'
    commandhelp = 'Filters PSMs by their confidence level. '

    def set_options(self):
        super().set_options()
        options = self.define_options(['confcol', 'conflvl', 'conftype',
                                       'unroll'], mzidtsv_options)
        self.options.update(options)

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.lowerbetter = self.conftype == 'lower'

    def get_psms(self):
        confkey = self.oldheader[int(self.confcol) - 1]
        self.header = self.oldheader[:]
        self.psms = prep.generate_psms(self.fn,
                                       self.oldheader,
                                       confkey,
                                       self.conflvl,
                                       self.lowerbetter)
