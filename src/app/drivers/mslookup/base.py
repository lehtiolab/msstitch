import os

from app.lookups import base as lookups
from app.drivers.base import BaseDriver
from app.drivers.options import mslookup_options


class LookupDriver(BaseDriver):
    outfile, outdir = None, None

    def __init__(self):
        super().__init__()
        self.parser_options = mslookup_options

    def set_options(self):
        super().set_options()
        del(self.options['-o'])
        del(self.options['-d'])
        self.options.update(self.define_options(['lookupfn'],
                                                mslookup_options))

    def initialize_lookup(self, outfile=None):
        if self.lookup is None:
            # FIXME MUST be a set or mzml lookup? here is place to assert
            # correct lookuptype!
            if outfile is None and self.outfile is None:
                self.outfile = os.path.join(self.outdir,
                                            'mslookup_db.sqlite')
                lookupfile = self.outfile
            elif outfile is not None:
                lookupfile = outfile
            elif self.outfile is not None:
                lookupfile = self.outfile
            self.lookup = lookups.create_new_lookup(lookupfile,
                                                    self.lookuptype)
        self.lookup.add_tables()

    def run(self):
        self.initialize_lookup()
        self.create_lookup()
