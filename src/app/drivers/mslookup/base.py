import os

from app.lookups import base as lookups
from app.drivers.base import BaseDriver
from app.drivers.options import mslookup_options


class LookupDriver(BaseDriver):
    def __init__(self):
        super().__init__()
        self.parser_options = mslookup_options

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['lookupfn'],
                                                mslookup_options))

    def initialize_lookup(self):
        if self.lookup is None:
            # FIXME MUST be a set or mzml lookup? here is place to assert
            # correct lookuptype!
            lookupfn = os.path.join(self.outdir,
                                    'mslookup_db.sqlite')
            self.lookup = lookups.create_new_lookup(lookupfn,
                                                    self.lookuptype)
        self.lookup.add_tables()

    def run(self):
        self.initialize_lookup()
        self.create_lookup()
