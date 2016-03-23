from app.actions.prottable import info as preparation
from app.drivers.prottable.base import ProttableAddData
from app.drivers.options import prottable_options
from app.writers import prottable as writers


class AddProteinInfoDriver(ProttableAddData):
    outsuffix = '_proteindata.txt'
    lookuptype = 'prottable'
    command = 'proteindata'
    commandhelp = ('Add protein data (description, coverage,# PSMs, etc.) '
                   'to a table with protein accessions')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['setname', 'genecentric',
                                                 'lookupfn'],
                                                prottable_options))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.headertypes = ['proteindata']
        self.poolnames = [kwargs.get('setname')]
        if self.genecentric:
            self.lookuptype = {'genes': 'genetable',
                               'assoc': 'associdtable'}[self.genecentric]

    def set_feature_generator(self):
        self.features = preparation.add_protein_data(self.in_proteins,
                                                     self.lookup,
                                                     self.headerfields,
                                                     self.genecentric,
                                                     self.setname
                                                     )

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_prottable_with_na(self.header, self.features, outfn)
