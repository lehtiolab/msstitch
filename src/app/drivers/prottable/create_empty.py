from app.actions.prottable import create_empty as preparation
from app.readers import tsv as reader
from app.drivers.prottable.base import ProttableDriver
from app.drivers.options import prottable_options


class CreateEmptyDriver(ProttableDriver):
    outsuffix = '_prottable.txt'
    command = 'emptytable'
    commandhelp = ('Create protein table from PSM table containing no '
                   'quant data, resulting in one column with (master) '
                   'proteins only. Use --protcol if input PSMs '
                   'are not of standard protein grouped variety.')

    def __init__(self):
        super().__init__()
        self.infiletype = 'TSV PSM table (MSGF+)'

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['proteincol']))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.headertypes = []

    def initialize_input(self):
        self.pepheader = reader.get_tsv_header(self.fn)
        self.in_psms = reader.generate_tsv_psms(self.fn, self.pepheader)

    def set_feature_generator(self):
        if self.proteincol is not None:
            self.get_column_header_for_number(['proteincol'], self.pepheader)
        self.features = preparation.generate_master_proteins(self.in_psms,
                                                             self.proteincol)
