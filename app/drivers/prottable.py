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
        self.precursorarea = False
        self.quantchannels = None
        self.oldheader = None

    def run(self):
        self.initialize_input()
        self.initialize_output()
        self.set_protein_generator()
        self.write()
        self.finish()

    def initialize_input(self):
        self.oldheader = reader.get_tsv_header(self.fn)
        self.in_proteins = reader.generate_tsv_proteins(self.fn, self.oldheader)

    def initialize_output(self):
        self.header = preparation.get_header(self.oldheader,
                                             self.quantchannels, self.protdata,
                                             self.precursorarea)

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_prottable(self.header, self.proteins, outfn)


class AddProteinInfoDriver(ProttableDriver):
    outsuffix = '_proteindata.txt'
    lookuptype = 'prottable'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def set_protein_generator(self):
        self.proteins = preparation.add_protein_data(self.in_proteins,
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
        self.isobaric = kwargs.get('isobaric', False)
        self.precursorarea = kwargs.get('precursor', False)

    def initialize_output(self):
        """Defines quantchannels from lookup table for header"""
        self.quantchannels = preparation.get_quantchannels(self.lookup)
        super().initialize_output()

    def initialize_input(self):
        """Not using input protein tables"""
        pass

    def set_protein_generator(self):
        """Generates proteins with quant from the lookup table"""
        self.proteins = preparation.build_quanted_proteintable(self.lookup,
                                                               self.header,
                                                               self.isobaric,
                                                               self.precursorarea)


class AddPrecursorAreaDriver(ProttableDriver):
    outsuffix = '_ms1q.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.precursorarea = True
        self.pepfile = kwargs.get('pepfile', False)

    def initialize_input(self):
        super().initialize_input()
        self.in_peptides = reader.generate_tsv_peptides(self.pepfile)

    def set_protein_generator(self):
        self.proteins = preparation.add_ms1_quant_from_top3_mzidtsv(
            self.in_proteins, self.in_peptides)
