import os

from app.actions.peptable import merge as preparation
from app.drivers.peptable.base import PeptableMergeDriver
from app.drivers.options import peptable_options


class BuildPeptideTableDriver(PeptableMergeDriver):
    outsuffix = ''
    lookuptype = 'peptidetable'
    command = 'build'
    commandhelp = ('Create peptide table from a lookup database. '
                   'E.g. when multiple peptide quant tables have '
                   'been read into the lookup and will be combined.'
                   'Does not need an -i.')

    def parse_input(self, **kwargs):
        """Build peptide table has no input file (though it has a lookup),
        which is why we set it to outfile name so the infile fetching
        and outfile creating wont error."""
        super().parse_input(**kwargs)
        self.fn = os.path.join(os.getcwd(), 'built_peptide_table.txt')
        if self.genecentric:
            self.lookuptype = 'peptidegenecentrictable'
        elif self.noncentric:
            self.genecentric = 'plain'
            self.lookuptype = 'peptidetableplain'

    def set_options(self):
        super().set_options()
        options = self.define_options(['lookupfn', 'genecentric', 'noncentric',
                                       'isobaric', 'precursor', 'fdr', 'pep',
                                       'mock_infn'], peptable_options)
        self.options.update(options)

    def set_feature_generator(self):
        """Generates proteins with quant from the lookup table"""
        self.features = preparation.build_peptidetable(self.lookup,
                                                       self.headerfields,
                                                       self.isobaric,
                                                       self.precursor,
                                                       self.fdr, self.pep,
                                                       self.genecentric)
