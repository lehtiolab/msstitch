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
                                                 'medianpsms', 'denompatterns',
                                                 'denomcols', 'minint'],
                                                mzidtsv_options))

    def get_psms(self):
        self.header = self.oldheader[:]
        if self.denomcols is not None:
            denomcols = [self.number_to_headerfield(col, self.oldheader)
                         for col in self.denomcols]
        elif self.denompatterns is not None:
            denomcolnrs = [tsv.get_columns_by_pattern(self.oldheader, pattern)
                           for pattern in self.denompatterns]
            denomcols = set([col for cols in denomcolnrs for col in cols])
        else:
            raise RuntimeError('Must define either denominator column numbers '
                               'or regex pattterns to find them')
        quantcols = tsv.get_columns_by_pattern(self.oldheader,
                                               self.quantcolpattern)
        if self.medianpsms is not None:
            medianheader = tsv.get_tsv_header(self.medianpsms)
        else:
            medianheader = False
        self.psms = prep.get_normalized_ratios(self.fn, self.oldheader,
                                               quantcols, denomcols,
                                               self.minint, self.medianpsms,
                                               medianheader)
