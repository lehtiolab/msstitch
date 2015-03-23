from app.actions.mslookup import proteingrouping as lookups
from app.drivers.mslookup import base
from app.readers import tsv as tsvreader


class ProteinGroupLookupDriver(base.LookupDriver):
    lookuptype = 'proteingroups'
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.spectracol = kwargs.get('speccol', None)

    def create_lookup(self):
        header = tsvreader.get_tsv_header(self.fn[0])
        confkey = header[int(self.confcol) - 1]
        specfncol = header[int(self.spectracol) - 1]
        lookups.create_protein_pep_lookup(self.fn,
                                          header,
                                          self.lookup,
                                          confkey,
                                          self.conflvl,
                                          self.lowerbetter,
                                          self.unroll,
                                          self.fasta,
                                          specfncol)
        lookups.build_proteingroup_db(self.fn, header,
                                      self.lookup,
                                      confkey, self.conflvl, self.lowerbetter,
                                      self.unroll, self.coverage)
