from app.actions.mslookup import proteinquant as lookups
from app.drivers.mslookup import base


class ProteinQuantLookupDriver(base.LookupDriver):
    """Creates lookup of protein tables that contain quant data"""
    lookuptype = 'prottable'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.quantcolpattern = kwargs.get('quantcolpattern', None)
        self.quantcolpattern = kwargs.get('psmnrcolpattern', None)
        self.proteincols = kwargs.get('protcol', None) - 1

    def create_lookup(self):
        lookups.create_proteinquant_lookup(self.fn, self.lookup,
                                           self.proteincols, 
                                           self.quantcolpattern,
                                           self.psmnrcolpattern)
