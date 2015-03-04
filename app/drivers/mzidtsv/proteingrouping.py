from app.preparation.mzidtsv import proteingrouping as prep
from app.drivers.mzidtsv import MzidTSVDriver
from app.readers import fasta
from app.lookups.sqlite.proteingroups import ProteinGroupProteinTableDB


class ProteinGroupDriver(MzidTSVDriver):
    outsuffix = '_protgroups.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.lookup = kwargs.get('lookup')

    def parse_fasta(self):
        if self.fasta:
            self.evidence_levels = fasta.has_evidence_levels(self.fasta)

    def get_psms(self):
        confkey = self.oldheader[int(self.confcol) - 1]
        protgroupdb = ProteinGroupProteinTableDB(self.lookup)
        self.header = prep.get_header_with_proteingroups(self.oldheader)
        self.psms = prep.generate_psms_with_proteingroups(self.fn,
                                                          self.oldheader,
                                                          self.header,
                                                          protgroupdb,
                                                          confkey,
                                                          self.conflvl,
                                                          self.lowerbetter,
                                                          self.unroll,
                                                          self.coverage,
                                                          self.evidence_levels)
