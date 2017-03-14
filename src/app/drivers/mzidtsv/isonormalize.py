from app.drivers.mzidtsv import MzidTSVDriver
from app.drivers.options import mzidtsv_options
from app.actions.mzidtsv import isonormalize as prep
from app.dataformats import prottable as prottabledata
from app.readers import tsv


class PSMIsoquantRatioDriver(MzidTSVDriver):
    outsuffix = '_ratio_isobaric.txt'
    command = 'isoratio'
    commandhelp = ('Produce isobaric PSM ratios from a '
                   'table containing raw intensities.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['quantcolpattern',
                                                 'denompatterns', 'denomcols',
                                                 'minint', 'targettable',
                                                 'proteincol', 'normalize',
                                                 'normalizeratios'],
                                                mzidtsv_options))

    def get_psms(self):
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
        self.get_column_header_for_number(['proteincol'], self.oldheader)
        nopsms = [prep.get_no_psms_field(qf) for qf in quantcols]
        if self.proteincol and self.targettable:
            targetheader = tsv.get_tsv_header(self.targettable)
            self.header = targetheader + quantcols + nopsms
        elif not self.proteincol and not self.targettable:
            self.header = (self.oldheader +
                           ['ratio_{}'.format(x) for x in quantcols])
        elif self.proteincol and not self.targettable:
            self.header = [prottabledata.HEADER_ACCESSION] + quantcols + nopsms
        self.psms = prep.get_isobaric_ratios(self.fn, self.oldheader,
                                             quantcols, denomcols, self.minint,
                                             self.targettable, self.proteincol,
                                             self.normalize,
                                             self.normalizeratios)


class PSMIsoquantNormalizeDriver(MzidTSVDriver):
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
