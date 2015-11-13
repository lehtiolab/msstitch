import os

from app.actions.peptable import merge as preparation
from app.drivers.peptable.base import PeptableMergeDriver


class BuildPeptideTableDriver(PeptableMergeDriver):
    outsuffix = ''
    lookuptype = 'peptidetable'

    def __init__(self, **kwargs):
        """Build peptide table has no input file (though it has a lookup),
        which is why we set it to outfile name so the infile fetching
        and outfile creating wont error."""
        kwargs['infile'] = os.path.join(os.getcwd(),
                                        'built_peptide_table.txt')
        self.genecentric = kwargs.get('genecentric', False)
        if self.genecentric:
            self.lookuptype = 'peptidegenecentrictable'
        super().__init__(**kwargs)
        self.isobaricquant = kwargs.get('isobaric', False)
        self.precursorquant = kwargs.get('precursor', False)
        self.fdr = kwargs.get('fdr', False)
        self.pep = kwargs.get('pep', False)

    def set_feature_generator(self):
        """Generates proteins with quant from the lookup table"""
        self.features = preparation.build_peptidetable(self.lookup,
                                                       self.header,
                                                       self.headerfields,
                                                       self.isobaricquant,
                                                       self.precursorquant,
                                                       self.fdr, self.pep,
                                                       self.genecentric)
