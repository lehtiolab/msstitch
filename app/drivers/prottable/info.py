from app.actions.prottable import info as preparation
from app.drivers.prottable.base import ProttableDriver


class AddProteinInfoDriver(ProttableDriver):
    outsuffix = '_proteindata.txt'
    lookuptype = 'prottable'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.protdata = True
        self.samplepool = kwargs.get('samplepool', False)

    def set_protein_generator(self):
        self.proteins = preparation.add_protein_data(self.in_proteins,
                                                     self.lookup,
                                                     self.samplepool)
