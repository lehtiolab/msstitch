from app.preparation.mzidtsv import proteingrouping as prep
from app.lookups import protein_peptide as lookups
from app.drivers.mzidtsv import MzidTSVDriver
from app.readers import fasta


class ProteinGroupDriver(MzidTSVDriver):
    outsuffix = '_protgroups.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.confcol = kwargs.get('confcol', None)
        self.conflvl = kwargs.get('conflvl', None)
        self.lowerbetter = kwargs.get('conftype', None) == 'lower'
        self.unroll = kwargs.get('unroll', False)
        self.evidence_levels = None
        self.fasta = kwargs.get('fasta', False)

    def parse_fasta(self):
        # FIXME maybe move this to preparation part.
        if self.fasta:
            self.evidence_levels = fasta.has_evidence_levels(self.fasta)

    def get_psms(self):
        if self.fasta:
            coverage = True
        confkey = self.oldheader[int(self.confcol) - 1]
        protgroupdb = lookups.create_protein_pep_lookup(self.fn,
                                                        self.oldheader,
                                                        confkey,
                                                        self.conflvl,
                                                        self.lowerbetter,
                                                        self.unroll,
                                                        self.fasta,
                                                        self.evidence_levels)
        self.header = prep.get_header_with_proteingroups(self.oldheader)
        self.psms = prep.generate_psms_with_proteingroups(self.fn,
                                                          self.oldheader,
                                                          self.header,
                                                          protgroupdb,
                                                          confkey,
                                                          self.conflvl,
                                                          self.lowerbetter,
                                                          self.unroll,
                                                          coverage,
                                                          self.evidence_levels)
