from app.preparation.mzidtsv import quant as prep
from app.readers import spectra as spectrareader
from app.readers import openms as openmsreader
from app.lookups import quant as lookups
from app.drivers.mzidtsv import MzidTSVDriver


class TSVQuantDriver(MzidTSVDriver):
    outsuffix = '_quant.tsv'

    def __init__(self, **kwargs):
        super(TSVQuantDriver, self).__init__(**kwargs)
        self.spectrafns = kwargs.get('spectra', None)
        self.quantfns = kwargs.get('quants', None)

    def get_psms(self):
        """Creates iterator to write to new tsv. Contains input tsv
        lines plus quant data for these."""
        quantdb = self.create_quantlookup()
        self.header, qheader = prep.get_full_and_quant_headers(
            self.oldheader, quantdb)
        self.psms = prep.generate_psms_quanted(quantdb, self.fn,
                                               qheader, self.oldheader)

    def create_quantlookup(self):
        """Creates sqlite file containing quantification data and
        spectra file/scan number data for looking up quant data.
        Returns sqlite file name.
        """
        fn_spectra = spectrareader.mzml_generator(self.spectrafns)
        consensus_quants = openmsreader.quant_generator(self.quantfns)
        return lookups.create_quant_lookup(fn_spectra, consensus_quants)
