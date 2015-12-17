from app.drivers.mslookup import base
from app.actions.mslookup import biosets as lookups


class BioSetLookupDriver(base.LookupDriver):
    outsuffix = '_setlookup.sqlite'
    lookuptype = 'biosets'
    command = 'biosets'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.setnames = kwargs.get('setnames')
        self.setnames = [x.replace('"', '') for x in self.setnames]

    def create_lookup(self):
        lookups.create_bioset_lookup(self.lookup, self.fn, self.setnames)
