import os

from app.actions.mzidtsv import merge_peptable as preparation
from app.drivers.prottable.base import ProttableMergeDriver


class BuildPeptideTableDriver(ProttableMergeDriver):
    outsuffix = ''
    lookuptype = 'peptidetable'

    def __init__(self, **kwargs):
        """Build peptide table has no input file (though it has a lookup),
        which is why we set it to outfile name so the infile fetching
        and outfile creating wont error."""
        kwargs['infile'] = os.path.join(os.getcwd(),
                                        'built_protein_table.txt')
        super().__init__(**kwargs)
        self.isobaricquant = kwargs.get('isobaric', False)
        self.precursorquant = kwargs.get('precursor', False)
        self.fdr = kwargs.get('fdr', False)
        self.pep = kwargs.get('pep', False)
        #self.proteindata = kwargs.get('proteindata', False)

    def set_protein_generator(self):
        """Generates proteins with quant from the lookup table"""
        self.psms = preparation.build_peptidetable(self.lookup,
                                                   self.header,
                                                   self.headerfields,
                                                   self.isobaricquant,
                                                   self.precursorquant,
                                                   self.fdr, self.pep)
                                                   #self.proteindata)
