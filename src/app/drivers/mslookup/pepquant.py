from app.actions.mslookup import proteinquant as lookups
from app.drivers.mslookup.proteinquant import ProteinQuantLookupDriver
from app.drivers.options import mslookup_options, peptable_options


class PeptideQuantLookupDriver(ProteinQuantLookupDriver):
    """Creates lookup of peptide tables that contain quant data"""
    lookuptype = 'peptidetable'
    command = 'peptides'

    def __init__(self):
        super().__init__()
        self.infiletype = 'peptide table'

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['peptidecol', 'setnames'],
                                                mslookup_options))
        self.options.update(self.define_options(['genecentric'],
                                                peptable_options))
        self.options.pop('--probcolpattern')
        self.options.pop('--protcol')

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        if self.genecentric:
            self.lookuptype = 'peptidegenecentrictable'

    def create_lookup(self):
        self.peptidecol = self.peptidecol - 1
        lookups.create_peptidequant_lookup(self.fn, self.lookup,
                                           self.setnames,
                                           self.peptidecol,
                                           self.precursorquantcolpattern,
                                           self.quantcolpattern,
                                           self.psmnrcolpattern,
                                           self.fdrcolpattern,
                                           self.pepcolpattern)
