from app.actions.prottable import create_from_psmtable as preparation
from app.readers import tsv as reader
from app.drivers.prottable.base import ProttableDriver


class CreateLabelfreeProteinDriver(ProttableDriver):
    outsuffix = '_labelfree.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.headertypes = ['precursorquant']

    def initialize_input(self):
        self.in_psms = reader.generate_tsv_peptides(self.fn)

    def set_protein_generator(self):
        self.proteins = preparation.create_protein_table_with_precursor_quant(
            self.in_psms, self.headerfields)
