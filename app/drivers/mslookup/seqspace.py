from app.drivers.mslookup import base
from app.actions.mslookup import searchspace as preparation
from app.drivers.options import mslookup_options


class SeqspaceLookupDriver(base.LookupDriver):
    """Creates an SQLite lookup DB from a FASTA file. Sequences are
    trypsinized and stored. It's possible to store sequences reversed
    for N-terminal falloff indexing, and it can be specified to cut
    tryptic before proline.
    """
    lookuptype = 'searchspace'
    command = 'seqspace'
    commandhelp = """Create a lookup DB from a FASTA file. Sequences are
    trypsinized and stored. It's possible to store sequences reversed
    for N-terminal falloff indexing, and it can be specified to cut
    tryptic before proline."""

    def set_options(self):
        super().set_options()
        self.options['--dbfile'].update({'required': False, 'default': None})
        self.options.update(self.define_options(['falloff', 'proline',
                                                 'trypsinize'],
                                                mslookup_options))

    def create_lookup(self):
        preparation.create_searchspace(self.lookup, self.fn, self.proline,
                                       self.falloff, self.trypsinize)
