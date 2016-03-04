from app.drivers.prottable.base import ProttableAddData
from app.actions.prottable import fdr as action
from app.readers import tsv as reader
from app.drivers.options import prottable_options


class ProttableFDRDriver(ProttableAddData):
    """Assigns FDR to protein table based on qvality output and protein table
    probabilities"""
    outsuffix = '_protfdr.txt'
    command = 'fdr'
    commandhelp = ('Add protein FDR to protein table by comparing '
                   'score (peptide q-value, protein probability, etc) '
                   'with qvality lookup table. Needs --scorecolpattern and '
                   'qvality output file specified with --qvality')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['qvalityout',
                                                 'scorecolpattern'],
                                                prottable_options))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.headertypes = ['proteinfdr', 'proteinpep']

    def initialize_input(self):
        super().initialize_input()
        self.scorecol = reader.get_cols_in_file(self.scorecolpattern,
                                                self.oldheader, True)

    def set_feature_generator(self):
        self.features = action.assign_protein_fdr(self.qvalityout,
                                                  self.in_proteins,
                                                  self.headerfields,
                                                  self.scorecol)
