import os

from app.actions.prottable import merge as preparation
from app.drivers.prottable.base import ProttableMergeDriver
from app.drivers.options import prottable_options


class BuildProteinTableDriver(ProttableMergeDriver):
    outsuffix = ''
    lookuptype = 'prottable'
    command = 'build'
    commandhelp = ('Create protein table from a lookup database. '
                   'E.g. when multiple protein quant tables have '
                   'been read into the lookup and will be combined.'
                   'Does not need an -i.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['lookupfn', 'genecentric',
                                                 'isobaric', 'precursor',
                                                 'probability', 'fdr', 'pep',
                                                 'mock_infn', 'mergecutoff'],
                                                prottable_options))

    def parse_input(self, **kwargs):
        """Build protein table has no input file (though it has a lookup),
        which is why we set it to outfile name so the infile fetching
        and outfile creating wont error."""
        super().parse_input(**kwargs)
        self.fn = os.path.join(os.getcwd(), 'built_protein_table.txt')
        if self.genecentric:
            self.lookuptype = {'genes': 'genetable',
                               'assoc': 'associdtable'}[self.genecentric]

    def set_feature_generator(self):
        """Generates proteins with quant from the lookup table"""
        self.features = preparation.build_proteintable(self.lookup,
                                                       self.headerfields,
                                                       self.mergecutoff,
                                                       self.isobaric,
                                                       self.precursor,
                                                       self.probability,
                                                       self.fdr, self.pep,
                                                       self.genecentric)
