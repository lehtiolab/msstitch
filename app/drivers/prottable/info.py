from app.actions.prottable import info as preparation
from app.drivers.prottable.base import ProttableAddData


class AddProteinInfoDriver(ProttableAddData):
    outsuffix = '_proteindata.txt'
    lookuptype = 'prottable'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.headertypes = ['proteindata']
        self.poolnames = [kwargs.get('setname')]
        self.genecentric = kwargs.get('genecentric', False)

    def set_feature_generator(self):
        self.features = preparation.add_protein_data(self.in_proteins,
                                                     self.lookup,
                                                     self.headerfields,
                                                     self.genecentric
                                                     )
