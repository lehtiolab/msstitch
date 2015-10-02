from app.drivers.mslookup import base
from app.actions.mslookup import searchspace as preparation


class SeqspaceLookupDriver(base.LookupDriver):
    """Creates an SQLite lookup DB from a FASTA file. Sequences are
    trypsinized and stored. It's possible to store sequences reversed
    for N-terminal falloff indexing, and it can be specified to cut
    tryptic before proline.
    """
    outsuffix = '_lookup.sqlite'
    lookuptype = 'searchspace'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.falloff = kwargs.get('falloff')
        self.proline = kwargs.get('proline')
        self.trypsinize = kwargs.get('notrypsin', True)

    def create_lookup(self):
        preparation.create_searchspace(self.lookup, self.fn[0], self.proline,
                                       self.falloff, self.trypsinize)
