from app.drivers.mslookup import base
from app.actions.mslookup import biosets as lookups
from app.drivers.options import mslookup_options


class BioSetLookupDriver(base.LookupDriver):
    outsuffix = '_setlookup.sqlite'
    lookuptype = 'biosets'
    command = 'biosets'
    commandhelp = ('Create SQLite lookup of mzML input files and '
                   'biological set names. Input files are passed to -i, '
                   'respective set names are passed to --setnames.')

    def __init__(self):
        super().__init__()
        self.infiletype = 'mzML spectra'

    def set_options(self):
        super().set_options()
        del(self.options['--dbfile'])
        self.options.update(self.define_options(['setnames'],
                                                mslookup_options))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.lookupfn = None
        self.setnames = [x.replace('"', '') for x in self.setnames]

    def create_lookup(self):
        lookups.create_bioset_lookup(self.lookup, self.fn, self.setnames)
