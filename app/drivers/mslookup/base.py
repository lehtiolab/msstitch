import os

from app.lookups import base as lookups
from app.drivers.base import BaseDriver


class LookupDriver(BaseDriver):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.lookupfn = kwargs.get('lookup', None)
        self.initialize_lookup()

    def initialize_lookup(self):
        if self.lookupfn is not None:
            self.lookup = lookups.get_lookup(self.lookupfn, self.lookuptype)
        else:
            # FIXME MUST be a set or mzml lookup? here is place to assert
            # correct lookuptype!
            self.lookupfn = os.path.join(self.outdir,
                                         'msstitcher_lookup.sqlite')
            self.lookup = lookups.create_new_lookup(self.lookupfn, self.lookuptype)
        self.lookup.add_tables()

    def run(self):
        self.create_lookup()
