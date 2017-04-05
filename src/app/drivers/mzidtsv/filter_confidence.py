from app.actions.mzidtsv import filter_confidence as prep
from app.drivers.mzidtsv import MzidTSVDriver
from app.drivers.options import mzidtsv_options
from app.readers import tsv as tsvreader


class ConfidenceFilterDriver(MzidTSVDriver):
    outsuffix = '_filtconf.txt'
    command = 'conffilt'
    commandhelp = 'Filters PSMs by their confidence level. '

    def set_options(self):
        super().set_options()
        options = self.define_options(['confcol', 'confpattern', 'conflvl',
                                       'conftype', 'unroll'], mzidtsv_options)
        self.options.update(options)

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.lowerbetter = self.conftype == 'lower'

    def get_psms(self):
        self.header = self.oldheader[:]
        if self.confpattern:
            confkey = tsvreader.get_cols_in_file(self.confpattern,
                                                 self.header, True)
        elif self.confcol:
            confkey = self.header[int(self.confcol) - 1]
        else:
            raise RuntimeError('Must define either --confcol or '
                               '--confcolpattern')
        self.psms = prep.generate_psms(self.fn,
                                       self.oldheader,
                                       confkey,
                                       self.conflvl,
                                       self.lowerbetter)
