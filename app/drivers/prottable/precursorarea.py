from app.readers import tsv as reader
from app.actions.prottable import precursorarea as preparation
from app.drivers.prottable.base import ProttableAddData
from app.drivers.options import prottable_options


class AddPrecursorAreaDriver(ProttableAddData):
    outsuffix = '_ms1q.txt'
    command = 'ms1quant'
    commandhelp = ('Add MS1 quantification data from a PSM table '
                   'containing precursor quant areas.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['psmfile', 'proteincol'],
                                                prottable_options))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.headertypes = ['precursorquant']

    def initialize_input(self):
        super().initialize_input()
        pepheader = reader.get_tsv_header(self.psmfile)
        self.get_column_header_for_number(['proteincol'], pepheader)
        self.in_peptides = reader.generate_tsv_peptides(self.psmfile)

    def set_feature_generator(self):
        self.features = preparation.add_ms1_quant_from_top3_mzidtsv(
            self.in_proteins, self.in_peptides, self.headerfields,
            self.proteincol)
