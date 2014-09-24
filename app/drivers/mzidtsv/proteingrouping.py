from app.preparation.mzidtsv import proteingrouping as prep
from app.lookups import protein_peptide as lookups
from app.drivers.mzidtsv import MzidTSVDriver


class ProteinGroupDriver(MzidTSVDriver):
    outsuffix = '_protgroups.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.confcol = kwargs.get('confcol', None)
        self.conflvl = kwargs.get('conflvl', None)
        self.lowerbetter = kwargs.get('conftype', None)
        self.unroll = kwargs.get('unroll', False)

    def get_psms(self):
        confkey = self.oldheader[self.conflvl]
        protgroupdb = lookups.create_protein_pep_lookup(self.fn,
                                                        self.oldheader,
                                                        confkey,
                                                        self.conflvl,
                                                        self.lowerbetter)
        self.header = prep.get_header_with_proteingroups(self.oldheader)
        self.psms = prep.generate_psms_with_proteingroups(self.fn,
                                                          self.oldheader,
                                                          protgroupdb,
                                                          confkey,
                                                          self.conflvl,
                                                          self.lowerbetter,
                                                          self.unroll)
