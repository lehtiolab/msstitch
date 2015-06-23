from app.drivers.prottable.base import ProttableDriver
from app.actions.prottable import probability as preparation
from app.readers import tsv as reader


class AddProteinProbability(ProttableDriver):
    outsuffix = '_protprob.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.pepfile = kwargs.get('pepfile', False)
        self.probability = True

    def initialize_input(self):
        super().initialize_input()
        self.in_peptides = reader.generate_tsv_peptides(self.pepfile)

    def set_protein_generator(self):
        self.proteins = preparation.add_nesvi_protein_probability(
            self.in_proteins, self.in_peptides)
