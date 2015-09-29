import os

from app.lookups import base as lookups
from app.drivers.base import BaseDriver


class LookupDriver(BaseDriver):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.conflvl = kwargs.get('conflvl', None)
        self.lowerbetter = kwargs.get('conftype', None) == 'lower'
        self.unroll = kwargs.get('unroll', False)
        self.fasta = kwargs.get('fasta', False)
        self.coverage = self.fasta is not False
        self.confcol = kwargs.get('confcol', False)
        self.proteincol = kwargs.get('protcol', False)

    def initialize_lookup(self):
        if self.lookup is None:
            # FIXME MUST be a set or mzml lookup? here is place to assert
            # correct lookuptype!
            lookupfn = os.path.join(self.outdir,
                                    'msstitcher_lookup.sqlite')
            self.lookup = lookups.create_new_lookup(lookupfn,
                                                    self.lookuptype)
        self.lookup.add_tables()

    def run(self):
        self.initialize_lookup()
        self.create_lookup()
