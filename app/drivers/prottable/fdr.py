from app.drivers.prottable.base import ProttableAddData
from app.actions.prottable import fdr as action


class ProttableFDRDriver(ProttableAddData):
    """Assigns FDR to protein table based on qvality output and protein table
    probabilities"""
    outsuffix = '_protfdr.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.qvalityfn = kwargs.get('qvalityfile')
        self.headertypes = ['proteinfdr', 'proteinpep']

    def set_protein_generator(self):
        self.proteins = action.assign_protein_fdr(self.qvalityfn,
                                                  self.in_proteins,
                                                  self.headerfields)
