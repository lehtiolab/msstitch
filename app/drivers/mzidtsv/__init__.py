from app.drivers import base
from app.readers import tsv as tsvreader
from app.writers import mzidtsv as writers
from app.drivers.options import mzidtsv_options


class MzidTSVDriver(base.BaseDriver):
    def __init__(self):
        super().__init__()
        self.parser_options = mzidtsv_options

    def run(self):
        self.oldheader = tsvreader.get_tsv_header(self.fn)
        self.get_psms()
        self.write()
        self.finish()

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_mzid_tsv(self.header, self.psms, outfn)
