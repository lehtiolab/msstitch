from app.readers import tsv as reader
from app.actions.prottable import precursorarea as preparation
from app.drivers.prottable.base import ProttableDriver


class AddPrecursorAreaDriver(ProttableDriver):
    outsuffix = '_ms1q.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.precursorarea = True
        self.pepfile = kwargs.get('psmfile', False)

    def initialize_input(self):
        super().initialize_input()
        self.in_peptides = reader.generate_tsv_peptides(self.pepfile)

    def set_protein_generator(self):
        self.proteins = preparation.add_ms1_quant_from_top3_mzidtsv(
            self.in_proteins, self.in_peptides)
