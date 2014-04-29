from app.drivers.basedrivers import PycolatorDriver
from app.lookups import sequences


class CreateLookup(PycolatorDriver):
    """Creates an SQLite lookup DB from a FASTA file. Sequences are
    trypsinized and stored. It's possible to store sequences reversed
    for N-terminal falloff indexing, and it can be specified to cut
    tryptic before proline.
    """
    outsuffix = '_lookup.sqlite'

    def __init__(self, **kwargs):
        super(PycolatorDriver, self).__init__(**kwargs)
        self.falloff = kwargs.get('falloff')
        self.proline = kwargs.get('proline')

    def run(self):
        self.outfn = self.create_outfilepath(self.fn, self.outsuffix)
        self.searchspace = sequences.create_searchspace(self.fn,
                                                        self.outfn,
                                                        self.proline,
                                                        self.falloff)
