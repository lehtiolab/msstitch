from app.writers import prottable as writers
from app.drivers.base import BaseDriver


class PepProttableDriver(BaseDriver):
    """Base class for prottable.py"""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.oldheader = False
        self.probability = False
        self.poolnames = False

    def run(self):
        self.initialize_input()
        self.create_header()
        self.set_feature_generator()
        self.write()
        self.finish()

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_prottable(self.header, self.features, outfn)

    def get_column_header_for_number(self, column_var_names, header=False):
        """This function subtracts 1 from inputted column number to comply
        with programmers counting (i.e. from 0, not from 1). Could possibly
        made more TSV general"""
        if not header:
            header = self.oldheader
        for col in column_var_names:
            value = getattr(self, col)
            if value is None:
                continue
            setattr(self, col, header[int(value) - 1])
