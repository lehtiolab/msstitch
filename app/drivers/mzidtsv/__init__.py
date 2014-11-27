from app.drivers import base
from app.readers import tsv as tsvreader
from app.writers import mzidtsv as writers


class MzidTSVDriver(base.BaseDriver):
    def run(self):
        self.oldheader = tsvreader.get_tsv_header(self.fn)
        self.parse_fasta()
        self.get_psms()
        self.write()
        self.finish()

    def parse_fasta(self):
        """Not implemented in all Drivers"""
        pass

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_mzid_tsv(self.header, self.psms, outfn)
