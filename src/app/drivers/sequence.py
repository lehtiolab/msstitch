from app.drivers import base
from app.actions import sequence
from app.actions.lookups import sequence as seqlup
from app.drivers.options import lookup_options, sequence_options
from Bio import SeqIO


class DecoySeqDriver(base.BaseDriver):
    outsuffix = '_decoy.fa'
    lookuptype = 'searchspace'
    command = 'makedecoy'
    commandhelp = """Create a decoy database from a FASTA file.
    tryp_rev reverses tryptic peptides
    prot_rev reverses full protein
    Both make sure no decoys are accidentally identical to a target sequence,
    unless --ignore-target-hits is passed"""

    def __init__(self):
        super().__init__()
        self.infiletype = 'FASTA'

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options([
            'fn', 'outfile', 'lookupfn', 'scramble', 'ignoretarget', 'trypsinize', 
            'miss_cleavage', 'minlength', 'max_shuffle'], sequence_options))
        self.options['lookupfn'].update({'required': False, 'default': None})

    def run(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        if self.lookup is None and not self.ignoretarget:
            self.initialize_lookup('decoychecker.sqlite')
            seqlup.create_searchspace(self.lookup, self.fn, self.minlength, reverse_seqs=False, miss_cleavage=self.miss_cleavage)
        decoyfa = sequence.create_decoy_fa(self.fn, self.scramble, self.lookup, self.trypsinize, self.miss_cleavage, self.minlength, self.max_shuffle)
        with open(outfn, 'w') as fp:
            SeqIO.write(decoyfa, fp, 'fasta')


class TrypsinizeDriver(base.BaseDriver):
    outsuffix = '_tryp.fa'
    command = 'trypsinize'
    commandhelp = """Trypsinize a FASTA file"""

    def __init__(self):
        super().__init__()
        self.infiletype = 'FASTA'

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options([
            'fn', 'outfile', 'miss_cleavage', 'minlength', 'proline'], sequence_options))

    def run(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        with open(self.fn) as fp, open(outfn, 'w') as wfp:
            seqs = SeqIO.parse(fp, 'fasta')
            decoyfa = sequence.create_trypsinized(seqs, self.proline, self.miss_cleavage, self.minlength)
            SeqIO.write(decoyfa, wfp, 'fasta')
