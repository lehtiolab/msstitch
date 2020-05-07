from app.actions.mzidtsv import proteingrouping as prep
from app.actions.mslookup import proteingrouping as lookupactions
from app.drivers.mzidtsv import MzidTSVDriver
from app.drivers.options import mzidtsv_options


class ProteinGroupDriver(MzidTSVDriver):
    outsuffix = '_protgroups.txt'
    lookuptype = 'proteingroups'
    command = 'proteingroup'
    commandhelp = ('Using a PSM table and corresponding lookup SQLite, '
            'this command adds calculated protein groups to it')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['lookupfn', 'unroll'],
                                                mzidtsv_options))

    def get_psms(self):
        self.lookup.add_tables()
        lookupactions.build_proteingroup_db(self.lookup)
        self.header = prep.get_header_with_proteingroups(self.oldheader)
        self.psms = prep.generate_psms_with_proteingroups(self.fn,
                                                          self.oldheader,
                                                          self.header,
                                                          self.lookup,
                                                          self.unroll)
