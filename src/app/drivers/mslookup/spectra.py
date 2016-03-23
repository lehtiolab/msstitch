from app.drivers.mslookup import base
from app.readers import spectra as spectrareader
from app.actions.mslookup import spectra as spectralookup
from app.actions.mslookup import biosets as biosetlookup
from app.drivers.options import mslookup_options


class SpectraLookupDriver(base.LookupDriver):
    lookuptype = 'spectra'
    command = 'spectra'
    commandhelp = ('Create lookup of spectra in mzML '
                   'format. Requires passing mzML files to -i, but '
                   'neither --spectra nor --dbfile. Biological set names '
                   'for each file should be specified using --setnames')

    def __init__(self):
        super().__init__()
        self.infiletype = 'mzML spectra'

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.spectrafns = self.fn
        if self.setnames is None:
            assert self.lookup is not None, ('Must specify lookup '
                                             'if setnames have not '
                                             'been provided')
        else:
            self.setnames = [x.replace('"', '') for x in self.setnames]

    def set_options(self):
        super().set_options()
        self.options['--dbfile'].update({'required': False, 'default': None})
        self.options.update(self.define_options(['multifiles', 'setnames'],
                                                mslookup_options))

    def create_lookup(self):
        biosetlookup.create_bioset_lookup(self.lookup, self.spectrafns,
                                          self.setnames)
        fn_spectra = spectrareader.mzmlfn_ms2_spectra_generator(
            self.spectrafns)
        spectralookup.create_spectra_lookup(self.lookup, fn_spectra)
