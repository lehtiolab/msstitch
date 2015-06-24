from app.readers import tsv as reader
from app.actions.prottable import old as preparation
from app.actions.prottable import headers
from app.writers import prottable as writers
from app.drivers.base import BaseDriver


class ProttableDriver(BaseDriver):
    """Base class for prottable.py"""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.protdata = kwargs.get('proteindata', False)
        self.precursorarea = False
        self.prottable_filenames = False
        self.quantchannels = False 
        self.oldheader = False
        self.probability = False
        self.samplepool = False

    def run(self):
        self.initialize_input()
        self.create_header()
        self.set_protein_generator()
        self.write()
        self.finish()

    def create_header(self):
        self.header = headers.get_header(self.oldheader,
                                         self.quantchannels, self.protdata,
                                         self.prottable_filenames,
                                         self.precursorarea, self.probability,
                                         self.samplepool)

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_prottable(self.header, self.proteins, outfn)


class ProttableAddData(ProttableDriver):
    def initialize_input(self):
        self.oldheader = reader.get_tsv_header(self.fn)
        self.in_proteins = reader.generate_tsv_proteins(self.fn, self.oldheader)


class ProttableMergeDriver(ProttableDriver):
    def initialize_input(self):
        self.quantchannels = preparation.get_quantchannels(self.lookup)
        self.prottable_filenames = preparation.get_precursorquant_headerfields(self.lookup)
