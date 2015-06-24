import os

from app.actions.prottable import merge as preparation
from app.drivers.prottable.base import ProttableMergeDriver


class BuildProteinTableDriver(ProttableMergeDriver):
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

    def set_protein_generator(self):
        """Generates proteins with quant from the lookup table"""
        self.proteins = preparation.build_proteintable(self.lookup,
                                                       self.header,
                                                       self.isobaric,
                                                       self.precursorarea,
                                                       self.probability)
