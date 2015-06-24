import os

from app.actions.prottable import merge as preparation
from app.drivers.prottable.base import ProttableDriver


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
        self.probability = True

    def initialize_output(self):
        """Defines quantchannels from lookup table for header"""
        self.quantchannels = preparation.get_quantchannels(self.lookup)
        self.prottable_filenames = preparation.get_precursorquant_headerfields(self.lookup)
        super().initialize_output()

    def initialize_input(self):
        """Not using input protein tables"""
        pass

    def set_protein_generator(self):
        """Generates proteins with quant from the lookup table"""
        self.proteins = preparation.build_proteintable(self.lookup,
                                                       self.header,
                                                       self.isobaric,
                                                       self.precursorarea,
                                                       self.probability)
