from app.actions.mzidtsv import spectra as actions
from app.drivers.mzidtsv import MzidTSVDriverLookup


class TSVSpectraDriver(MzidTSVDriverLookup):
    lookuptype = 'spectra'
    outsuffix = '_spectradata.tsv'

    def get_psms(self):
        """Creates iterator to write to new tsv. Contains input tsv
        lines plus quant data for these."""
        self.header = actions.create_header(self.oldheader, self.spec_column)
        self.psms = actions.generate_psms_spectradata(self.lookup, self.fn,
                                               self.oldheader, self.spec_column)
