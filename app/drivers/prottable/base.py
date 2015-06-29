from app.readers import tsv as reader
from app.actions.prottable import headers
from app.writers import prottable as writers
from app.drivers.base import BaseDriver


class ProttableDriver(BaseDriver):
    """Base class for prottable.py"""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.oldheader = False
        self.probability = False
        self.poolnames = False

    def run(self):
        self.initialize_input()
        self.create_header()
        self.set_protein_generator()
        self.write()
        self.finish()

    def create_header(self):
        self.headerfields = headers.get_headerfields(self.headertypes, self.lookup,
                                                     self.poolnames)
        self.header = headers.generate_header(self.headerfields, self.oldheader)

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_prottable(self.header, self.proteins, outfn)


class ProttableAddData(ProttableDriver):
    def initialize_input(self):
        self.oldheader = reader.get_tsv_header(self.fn)
        self.in_proteins = reader.generate_tsv_proteins(self.fn, self.oldheader)


class ProttableMergeDriver(ProttableDriver):
    def initialize_input(self):
        self.headertypes = []
        for inflag, htype in zip([self.proteindata, self.probability, self.precursorquant, self.isobaricquant], ['proteindata', 'probability', 'precursorquant', 'isoquant']):
            if inflag:
                self.headertypes.append(htype)
        self.poolnames = [x[0] for x in self.lookup.get_all_poolnames()]
