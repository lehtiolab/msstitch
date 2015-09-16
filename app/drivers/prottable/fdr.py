from app.drivers.prottable.base import ProttableAddData
from app.actions.prottable import fdr as action


class ProttableFDRDriver(ProttableAddData):
    """Assigns FDR to protein table based on qvality output and protein table
    probabilities"""
    outsuffix = '_protfdr.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.qvalityfn = kwargs.get('qvalityfile')
        self.scorecol = kwargs.get('scorecol', None)
        self.headertypes = ['proteinfdr', 'proteinpep']

    def initialize_input(self):
        super().initialize_input()
        self.get_column_header_for_number(['scorecol'])

    def set_feature_generator(self):
        self.features = action.assign_protein_fdr(self.qvalityfn,
                                                  self.in_proteins,
                                                  self.headerfields,
                                                  self.scorecol)
