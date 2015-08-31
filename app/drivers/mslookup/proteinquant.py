from app.actions.mslookup import proteinquant as lookups
from app.drivers.mslookup import base


class ProteinQuantLookupDriver(base.LookupDriver):
    """Creates lookup of protein tables that contain quant data"""
    lookuptype = 'prottable'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.poolnames = [x.replace('"', '') for x in kwargs.get('setnames')]
        # FIXME need check to see same poolnames correlate with self.fn len
        self.quantcolpattern = kwargs.get('quantcolpattern', None)
        self.psmnrcolpattern = kwargs.get('psmnrcolpattern', None)
        self.precursorquantcolpattern = kwargs.get('precursorquantcolpattern',
                                                   None)
        self.proteincols = kwargs.get('protcol', None) - 1
        self.probcolpattern = kwargs.get('probcolpattern', None)

    def create_lookup(self):
        lookups.create_proteinquant_lookup(self.fn, self.lookup,
                                           self.poolnames,
                                           self.proteincols,
                                           self.precursorquantcolpattern,
                                           self.quantcolpattern,
                                           self.psmnrcolpattern,
                                           self.probcolpattern)
