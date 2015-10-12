from app.actions.prottable import create_empty as preparation
from app.readers import tsv as reader
from app.drivers.prottable.base import ProttableDriver


class CreateEmptyDriver(ProttableDriver):
    outsuffix = '.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.headertypes = []
        self.protcol = kwargs.get('protcol')

    def initialize_input(self):
        self.pepheader = reader.get_tsv_header(self.fn)
        self.in_psms = reader.generate_tsv_psms(self.fn, self.pepheader)

    def set_feature_generator(self):
        if self.protcol is not None:
            self.get_column_header_for_number(['protcol'], self.pepheader)
        self.features = preparation.generate_master_proteins(self.in_psms, self.protcol)
