from app.actions.prottable import create_empty as preparation
from app.readers import tsv as reader
from app.drivers.prottable.base import ProttableDriver


class CreateEmptyDriver(ProttableDriver):
    outsuffix = '.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.headertypes = []

    def initialize_input(self):
        self.in_psms = reader.generate_tsv_peptides(self.fn)

    def set_protein_generator(self):
        self.proteins = preparation.generate_master_proteins(self.in_psms)
