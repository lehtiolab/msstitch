from app.writers import prottable as writers
from app.drivers.base import BaseDriver


class PepProttableDriver(BaseDriver):
    """Base class for prottable.py"""
    def __init__(self):
        super().__init__()
        self.oldheader = False
        self.probability = False
        self.poolnames = False
        self.group_by_field = False

    def run(self):
        self.initialize_input()
        self.create_header()
        self.set_feature_generator()
        self.write()
        self.finish()

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_prottable(self.header, self.features, outfn)
