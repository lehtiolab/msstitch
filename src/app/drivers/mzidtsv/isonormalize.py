from app.drivers.mzidtsv import MzidTSVDriver
from app.drivers.options import mzidtsv_options
from app.actions.mzidtsv import isonormalize as prep
from app.readers import tsv


class MzidTSVIsoquantNormalizeDriver(MzidTSVDriver):
    outsuffix = '_normalized_isobaric.txt'
    command = 'isonormalize'
    commandhelp = ('Produce median normalized isobaric PSM ratios from a '
                   'table containing raw intensities.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['quantcolpattern',
                                                 'medianpsms',
                                                 'denomcols', 'minint'],
                                                mzidtsv_options))

    def get_psms(self):
        self.header = self.oldheader[:]
        denomcols = [self.number_to_headerfield(col, self.oldheader)
                     for col in self.denomcols]
        quantcols = tsv.get_columns_by_pattern(self.oldheader,
                                               self.quantcolpattern)
        self.psms = prep.get_normalized_ratios(self.fn, self.oldheader,
                                               quantcols, denomcols,
                                               self.minint, self.medianpsms)
