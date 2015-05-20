import os

from app.drivers.base import BaseDriver
from app.writers import prottable as writers
from app.readers import tsv as reader
from app.actions import prottable as preparation


class ProttableDriver(BaseDriver):
    """Base class for prottable.py"""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.protdata = kwargs.get('proteindata', False)
        self.quantchannels = None
        self.oldheader = None

    def run(self):
        self.initialize_output()
        self.set_protein_generator()
        self.write()
        self.finish()

    def initialize_output(self):
        self.header = preparation.get_header(self.oldheader, self.quantchannels, self.protdata)

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_prottable(self.header, self.proteins, outfn)


class AddProteinInfoDriver(ProttableDriver):
    outsuffix = '_proteindata.txt'
    lookuptype = 'prottable'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.oldheader = reader.get_tsv_header(self.fn)

    def set_protein_generator(self):
        proteins = reader.generate_tsv_proteins(self.fn, self.oldheader)
        self.proteins = preparation.add_protein_data(proteins,
                                                     self.lookup)


class BuildProteinTableDriver(ProttableDriver):
    outsuffix = ''
    lookuptype = 'prottable'

    def __init__(self, **kwargs):
        """Build protein table has no input file (though it has a lookup),
        which is why we set it to outfile name so the infile fetching
        and outfile creating wont error."""
        kwargs['infile'] = os.path.join(os.getcwd(),
                                        'built_protein_table.txt')
        super().__init__(**kwargs)

    def initialize_output(self):
        """Defines quantchannels from lookup table for header"""
        self.quantchannels = preparation.get_quantchannels(self.lookup)
        super().initialize_output()

    def set_protein_generator(self):
        """Generates proteins with quant from the lookup table"""
        self.proteins = preparation.build_quanted_proteintable(self.lookup,
                                                               self.header)
