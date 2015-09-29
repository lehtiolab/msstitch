from app.actions.mslookup import proteingrouping as lookups
from app.drivers.mslookup import base
from app.readers import tsv as tsvreader


class ProteinGroupLookupDriver(base.LookupDriver):
    lookuptype = 'proteingroups'

    def create_lookup(self):
        header = tsvreader.get_tsv_header(self.fn[0])
        self.get_column_header_for_number(['confcol', 'proteincol'], header)
        allpsms = lookups.create_protein_pep_lookup(self.fn,
                                                    header,
                                                    self.lookup,
                                                    self.confcol,
                                                    self.conflvl,
                                                    self.lowerbetter,
                                                    self.fasta,
                                                    self.proteinfield)
        lookups.build_proteingroup_db(self.lookup, allpsms, self.coverage)
