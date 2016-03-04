from app.drivers import base
from app.readers import tsv as tsvreader
from app.writers import mzidtsv as writers
from app.drivers.options import mzidtsv_options


class MzidTSVDriver(base.BaseDriver):
    def __init__(self):
        super().__init__()
        self.infiletype = 'TSV PSM table (MSGF+)'
        self.parser_options = mzidtsv_options

    def run(self):
        if type(self.fn) == list:
            self.first_infile = self.fn[0]
        else:
            self.first_infile = self.fn
        self.oldheader = tsvreader.get_tsv_header(self.first_infile)
        self.get_psms()
        self.write()
        self.finish()

    def write(self):
        outfn = self.create_outfilepath(self.first_infile, self.outsuffix)
        writers.write_mzid_tsv(self.header, self.psms, outfn)
