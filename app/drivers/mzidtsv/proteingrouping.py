import shutil

from app.preparation.mzidtsv import proteingrouping as prep
from app.lookups import protein_peptide as lookups
from app.drivers.mzidtsv import MzidTSVDriver
from app.readers import fasta


class ProteinGroupDriver(MzidTSVDriver):
    outsuffix = '_protgroups.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.confcol = kwargs.get('confcol', None)
        self.conflvl = kwargs.get('conflvl', None)
        self.lowerbetter = kwargs.get('conftype', None) == 'lower'
        self.unroll = kwargs.get('unroll', False)
        self.evidence_levels = None
        self.fasta = kwargs.get('fasta', False)
        self.coverage = self.fasta is not False

    def parse_fasta(self):
        if self.fasta:
            self.evidence_levels = fasta.has_evidence_levels(self.fasta)

    def create_protein_pep_lookup(self, confkey):
        return lookups.create_protein_pep_lookup(self.fn,
                                                 self.workdir,
                                                 self.oldheader,
                                                 confkey,
                                                 self.conflvl,
                                                 self.lowerbetter,
                                                 self.unroll,
                                                 self.fasta,
                                                 self.evidence_levels)

    def get_psms(self):
        confkey = self.oldheader[int(self.confcol) - 1]
        protgroupdb = self.create_protein_pep_lookup(confkey)
        self.header = prep.get_header_with_proteingroups(self.oldheader)
        self.psms = prep.generate_psms_with_proteingroups(self.fn,
                                                          self.oldheader,
                                                          self.header,
                                                          protgroupdb,
                                                          confkey,
                                                          self.conflvl,
                                                          self.lowerbetter,
                                                          self.unroll,
                                                          self.coverage,
                                                          self.evidence_levels)


class ProteinGroupLookupDriver(ProteinGroupDriver):
    outsuffix = '_protgrouplookup.txt'

    def get_psms(self):
        # FIXME Should really be called process_psms, but have to change in
        # all mzidtsv drivers
        confkey = self.oldheader[int(self.confcol) - 1]
        self.prot_pep_lookup = self.create_protein_pep_lookup(confkey)
        prep.build_proteingroup_db(self.fn, self.oldheader,
                                   self.prot_pep_lookup,
                                   confkey, self.conflvl, self.lowerbetter,
                                   self.unroll, self.coverage)

    def write(self):
        """Moves outfile from workdir to destination"""
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        shutil.move(self.prot_pep_lookup.get_fn(), outfn)
