from app.drivers.prottable.base import ProttableAddData
from app.actions.prottable import probability as preparation
from app.readers import tsv as reader


class AddProteinProbability(ProttableAddData):
    outsuffix = '_protprob.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.pepfile = kwargs.get('pepfile', False)
        self.headertypes = ['probability']

    def initialize_input(self):
        super().initialize_input()
        self.in_peptides = reader.generate_tsv_peptides(self.pepfile)

    def set_feature_generator(self):
        self.features = preparation.add_nesvi_protein_probability(
            self.in_proteins, self.in_peptides, self.headerfields)
