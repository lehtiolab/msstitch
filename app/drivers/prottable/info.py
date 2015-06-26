from app.actions.prottable import info as preparation
from app.drivers.prottable.base import ProttableAddData


class AddProteinInfoDriver(ProttableAddData):
    outsuffix = '_proteindata.txt'
    lookuptype = 'prottable'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.headertypes = ['proteindata']
        self.poolnames = [kwargs.get('setname')]

    def set_protein_generator(self):
        self.proteins = preparation.add_protein_data(self.in_proteins,
                                                     self.lookup,
                                                     self.headerfields,
                                                     )
