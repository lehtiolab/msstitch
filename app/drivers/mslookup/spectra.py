from app.drivers.mslookup import base
from app.readers import spectra as spectrareader
from app.actions.mslookup import spectra as lookups


class SpectraLookupDriver(base.LookupDriver):
    outsuffix = '_spectralookup.sqlite'
    lookuptype = 'spectra'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.spectrafns = self.fn
        self.setnames = kwargs.get('setnames', None)
        if self.setnames is None:
            assert self.lookup is not None, ('Must specify lookup '
                                             'if setnames have not '
                                             'been provided')

    def create_lookup(self):
        lookups.create_bioset_lookup(self.lookup, self.spectrafns,
                                     self.setnames)
        fn_spectra = spectrareader.mzmlfn_spectra_generator(self.spectrafns)
        lookups.create_spectra_lookup(self.lookup, fn_spectra)
