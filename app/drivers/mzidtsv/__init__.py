from app.drivers import base
from app.readers import tsv as tsvreader
from app.writers import mzidtsv as writers
from app.lookups import base as lookups


class MzidTSVDriver(base.BaseDriver):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.confcol = kwargs.get('confcol', None)
        self.conflvl = kwargs.get('conflvl', None)
        self.lowerbetter = kwargs.get('conftype', None) == 'lower'
        self.unroll = kwargs.get('unroll', False)
        self.evidence_levels = None
        self.fasta = kwargs.get('fasta', False)
        self.coverage = self.fasta is not False
        self.spec_column = kwargs.get('speccol', None)

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
