from app.actions.mslookup import proteinquant as lookups
from app.drivers.mslookup import base


class ProteinQuantLookupDriver(base.LookupDriver):
    """Creates lookup of protein tables that contain quant data"""
    lookuptype = 'prottable'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.poolnames = [x.replace('"', '') for x in kwargs.get('setnames')]
        # FIXME need check to see same poolnames correlate with self.fn len
        self.quantcolpattern = kwargs.get('isobquantcolpattern', None)
        self.psmnrcolpattern = kwargs.get('psmnrcolpattern', None)
        self.precursorquantcolpattern = kwargs.get('precursorquantcolpattern',
                                                   None)
        self.proteincol = kwargs.get('protcol', None)
        self.probcolpattern = kwargs.get('probcolpattern', None)
        self.fdrcolpattern = kwargs.get('fdrcolpattern', None)
        self.pepcolpattern = kwargs.get('pepcolpattern', None)

    def create_lookup(self):
        self.proteincol = self.proteincol - 1
        lookups.create_proteinquant_lookup(self.fn, self.lookup,
                                           self.poolnames,
                                           self.proteincol,
                                           self.precursorquantcolpattern,
                                           self.quantcolpattern,
                                           self.psmnrcolpattern,
                                           self.probcolpattern,
                                           self.fdrcolpattern,
                                           self.pepcolpattern)
