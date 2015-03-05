from app.actions.mslookup import proteingrouping as lookups
from app.drivers.mslookup import base


class ProteinGroupLookupDriver(base.LookupDriver):

    def create_lookup(self):
        confkey = self.oldheader[int(self.confcol) - 1]
        lookups.create_protein_pep_lookup(self.fn,
                                          self.oldheader,
                                          confkey,
                                          self.conflvl,
                                          self.lowerbetter,
                                          self.unroll,
                                          self.fasta,
                                          self.evidence_levels)
        lookups.build_proteingroup_db(self.fn, self.oldheader,
                                      self.lookup,
                                      confkey, self.conflvl, self.lowerbetter,
                                      self.unroll, self.coverage)
