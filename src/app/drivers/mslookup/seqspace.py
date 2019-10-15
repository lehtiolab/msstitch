from app.drivers.mslookup import base
from app.actions.mslookup import searchspace as preparation
from app.drivers.options import mslookup_options, sequence_options
from Bio import SeqIO


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

    def __init__(self):
        super().__init__()
        self.infiletype = 'FASTA'

    def set_options(self):
        super().set_options()
        self.options['--dbfile'].update({'required': False, 'default': None})
        self.options.update(self.define_options(['falloff', 'proline'], mslookup_options))
        self.options.update(self.define_options(['trypsinize'], sequence_options))

    def create_lookup(self):
        preparation.create_searchspace(self.lookup, self.fn, self.proline,
                                       self.falloff, self.trypsinize)


class WholeProteinSeqspaceLookupDriver(base.LookupDriver):
    lookuptype = 'searchspace'
    command = 'protspace'
    commandhelp = """Create a full-length protein lookup DB from a FASTA file.
    Protein sequences are stored as short peptides starting from every one of
    the proteins' amino acids."""

    def __init__(self):
        super().__init__()
        self.infiletype = 'FASTA'

    def set_options(self):
        super().set_options()
        self.options['--dbfile'].update({'required': False, 'default': None})
        self.options.update(self.define_options(['minlength'],
                                                mslookup_options))

    def create_lookup(self):
        preparation.create_searchspace_wholeproteins(self.lookup, self.fn,
                                                     self.minlength)



class DecoySeqDriver(base.LookupDriver):
    outsuffix = '_decoy.fa'
    lookuptype = 'searchspace'
    command = 'makedecoy'
    commandhelp = """Create a decoy database from a FASTA file.
    tryp_rev reverses tryptic peptides
    prot_rev reverses full protein
    Both make sure no decoys are accidentally identical to a target sequence,
    unless --ignore-target-hits is passed"""
    # FIXME doesnt really fit in mslookup category
    # maybe we should dump the categories in next version

    def __init__(self):
        super().__init__()
        self.infiletype = 'FASTA'

    def set_options(self):
        super().set_options()
        self.options['--dbfile'].update({'required': False, 'default': None})
        self.options.update(self.define_options(
            ['fn', 'outfile', 'scramble', 'ignoretarget', 'trypsinize'], sequence_options))

    def run(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        if self.lookup is None and not self.ignoretarget:
            self.initialize_lookup('decoychecker.sqlite')
            preparation.create_searchspace(self.lookup, self.fn, reverse_seqs=False, fully_tryptic=True)
        decoyfa = preparation.create_decoy_fa(self.fn, self.scramble, self.lookup, self.trypsinize)
        with open(outfn, 'w') as fp:
            SeqIO.write(decoyfa, fp, 'fasta')
