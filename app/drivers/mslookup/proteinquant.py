from app.actions.mslookup import proteinquant as lookups
from app.drivers.mslookup import base


class ProteinQuantLookupDriver(base.LookupDriver):
    """Creates lookup of protein tables that contain quant data"""
    lookuptype = 'proteinquant'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.quantcols = kwargs.get('quantcols', None)

    def create_lookup(self):
        lookups.create_proteinquant_lookup(self.fn, self.lookup,
                                           self.quantcols)
