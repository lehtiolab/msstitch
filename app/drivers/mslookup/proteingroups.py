from app.lookups import protein_peptide as lookups
from app.drivers.mslookup import base
from app.preparation.mzidtsv import proteingrouping as prep
# FIXME move from mzidtsv to mslookup


class ProteinGroupLookupDriver(base.LookupDriver):
    outsuffix = '_protgrouplookup.txt'

    def create_protein_pep_lookup(self, confkey):
        return lookups.create_protein_pep_lookup(self.fn,
                                                 self.oldheader,
                                                 confkey,
                                                 self.conflvl,
                                                 self.lowerbetter,
                                                 self.unroll,
                                                 self.fasta,
                                                 self.evidence_levels)

    def get_psms(self):
        # FIXME Should really be called process_psms, but have to change in
        # all mzidtsv drivers
        confkey = self.oldheader[int(self.confcol) - 1]
        self.lookup = self.create_protein_pep_lookup(confkey)
        prep.build_proteingroup_db(self.fn, self.oldheader,
                                   self.lookup,
                                   confkey, self.conflvl, self.lowerbetter,
                                   self.unroll, self.coverage)
