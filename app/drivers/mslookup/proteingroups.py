from app.actions.mslookup import proteingrouping as lookups
from app.drivers.mslookup import base
from app.readers import tsv as tsvreader


class ProteinGroupLookupDriver(base.LookupDriver):
    lookuptype = 'proteingroups'

    def create_lookup(self):
        header = tsvreader.get_tsv_header(self.fn[0])
        self.get_column_header_for_number(['proteincol'], header)
        lookups.build_proteingroup_db(self.lookup)
