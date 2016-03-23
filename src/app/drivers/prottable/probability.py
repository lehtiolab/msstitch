from app.drivers.prottable.base import ProttableAddData
from app.actions.prottable import probability as preparation
from app.readers import tsv as reader
from app.drivers.options import prottable_options


class AddProteinProbability(ProttableAddData):
    outsuffix = '_protprob.txt'
    command = 'probability'
    commandhelp = ('Add protein probabilities from peptide table posterior '
                   'error probabilities (PEP). Probabilities are calculated '
                   'as in Nesvizhskii et al. (2003) Anal.Chem., eq 3.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['pepfile', 'proteincol'],
                                                prottable_options))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.headertypes = ['probability']

    def initialize_input(self):
        super().initialize_input()
        self.in_peptides = reader.generate_tsv_peptides(self.pepfile)
        pepheader = reader.get_tsv_header(self.pepfile)
        self.get_column_header_for_number(['proteincol'], pepheader)

    def set_feature_generator(self):
        self.features = preparation.add_nesvi_protein_probability(
            self.in_proteins, self.in_peptides, self.headerfields,
            self.proteincol)
