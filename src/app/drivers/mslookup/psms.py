from app.drivers.mslookup import base
from app.readers import tsv as tsvreader
from app.actions.mslookup import psms as lookup
from app.drivers.options import mslookup_options


class PSMLookupDriver(base.LookupDriver):
    lookuptype = 'psm'
    command = 'psms'
    commandhelp = ('Loads PSM table into lookup. Important for '
                   'several steps later on, such as protein grouping and '
                   'PSM quantitation. PSM TSV table passed to -i, '
                   'With flag --unroll, --spectracol, and specify a FASTA '
                   'file with --fasta, a Biomart map file with --map. When '
                   'using --map for a decoy lookup, use --decoy.')

    def __init__(self):
        super().__init__()
        self.infiletype = 'TSV PSM table (MSGF+)'

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['spectracol', 'mapfn',
                                                 'decoy', 'fasta', 'unroll',
                                                 'genefield', 'fastadelim'],
                                                mslookup_options))

    def create_lookup(self):
        header = tsvreader.get_tsv_header(self.fn)
        specfncol = header[int(self.spectracol) - 1]
        fastadelim, genefield = self.get_fastadelim_genefield(self.fastadelim,
                                                              self.genefield)
        lookup.create_psm_lookup(self.fn, self.fasta, self.mapfn, header,
                                 self.lookup, self.unroll, specfncol,
                                 self.decoy, fastadelim, genefield)
