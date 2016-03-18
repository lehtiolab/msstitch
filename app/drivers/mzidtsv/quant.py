from app.actions.mzidtsv import quant as prep
from app.drivers.mzidtsv import MzidTSVDriver
from app.drivers.options import mzidtsv_options


class TSVQuantDriver(MzidTSVDriver):
    lookuptype = 'quant'
    outsuffix = '_quant.tsv'
    command = 'quant'
    commandhelp = ('Add quantitative data to a tab separated file with'
                   'PSMs. Quant data is fetched from an SQLite lookup '
                   'Specify --isobaric and/or --precursor.')

    def set_options(self):
        super().set_options()
        options = self.define_options(['lookupfn', 'precursor', 'isobaric'],
                                      mzidtsv_options)
        self.options.update(options)

    def get_psms(self):
        """Creates iterator to write to new tsv. Contains input tsv
        lines plus quant data for these."""
        self.header, isob_header = prep.get_full_and_isobaric_headers(
            self.oldheader, self.lookup, self.isobaric, self.precursor)
        self.psms = prep.generate_psms_quanted(self.lookup, self.fn,
                                               isob_header, self.oldheader,
                                               self.isobaric, self.precursor)
