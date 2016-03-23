from app.actions.mzidtsv import proteingrouping as prep
from app.drivers.mzidtsv import MzidTSVDriver
from app.drivers.options import mzidtsv_options


class ProteinGroupDriver(MzidTSVDriver):
    outsuffix = '_protgroups.txt'
    lookuptype = 'proteingroups'
    command = 'proteingroup'
    commandhelp = ('Takes lookup SQLite result, uses it to output '
                   'mzidtsv file with protein groups')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['lookupfn', 'unroll'],
                                                mzidtsv_options))

    def get_psms(self):
        self.header = prep.get_header_with_proteingroups(self.oldheader)
        self.psms = prep.generate_psms_with_proteingroups(self.fn,
                                                          self.oldheader,
                                                          self.header,
                                                          self.lookup,
                                                          self.unroll)
