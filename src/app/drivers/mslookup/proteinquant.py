from app.actions.mslookup import proteinquant as lookups
from app.drivers.mslookup import base
from app.drivers.options import mslookup_options


class ProteinQuantLookupDriver(base.LookupDriver):
    """Creates lookup of protein tables that contain quant data"""
    lookuptype = 'prottable'
    command = 'proteins'
    commandhelp = ('Creates lookup of protein quantification '
                   'data in tab separated format. Header should include '
                   'quantification channel names, and if possible the '
                   'number of peptides quantified for each protein in the '
                   'respective channels. Can be used with --genecentric '
                   'Lookup should already include proteins.')

    def __init__(self):
        super().__init__()
        self.infiletype = 'protein table'

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['genecentric', 'setnames',
                                                 'quantcolpattern',
                                                 'psmnrcolpattern',
                                                 'precursorquantcolpattern',
                                                 'probcolpattern',
                                                 'fdrcolpattern',
                                                 'pepcolpattern',
                                                 'multifiles', 'proteincol'],
                                                mslookup_options))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        if self.genecentric in ['genes', 'assoc']:
            # cannot simply do if genecentric: because sometimes it is
            # True e.g. from a child class and then a KeyError arises below
            self.lookuptype = {'genes': 'genetable',
                               'assoc': 'associdtable'}[self.genecentric]
        self.setnames = [x.replace('"', '') for x in self.setnames]
        # FIXME need check to see same poolnames correlate with self.fn len

    def create_lookup(self):
        self.proteincol = self.proteincol - 1
        lookups.create_proteinquant_lookup(self.fn, self.lookup,
                                           self.setnames,
                                           self.proteincol,
                                           self.precursorquantcolpattern,
                                           self.quantcolpattern,
                                           self.psmnrcolpattern,
                                           self.probcolpattern,
                                           self.fdrcolpattern,
                                           self.pepcolpattern)
