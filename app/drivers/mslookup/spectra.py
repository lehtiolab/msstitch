from app.drivers.mslookup import base
from app.readers import spectra as spectrareader
from app.actions.mslookup import spectra as lookups


class SpectraLookupDriver(base.LookupDriver):
    outsuffix = '_spectralookup.sqlite'
    lookuptype = 'spectra'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.spectrafns = self.fn

    def create_lookup(self):
        fn_spectra = spectrareader.mzmlfn_spectra_generator(self.spectrafns)
        lookups.create_spectra_lookup(self.lookup, fn_spectra)
