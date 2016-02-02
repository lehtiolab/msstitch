from app.actions.mzidtsv import prot2gene as actions
from app.drivers.mzidtsv import MzidTSVDriver
from app.drivers.options import mzidtsv_options


class TSVGeneFromProteinDriver(MzidTSVDriver):
    outsuffix = '_genes.txt'
    lookuptype = 'psm'
    command = 'genes'
    commandhelp = ('Add column to mzidtsv with gene names or symbols, '
                   'which are stored in a lookup specified with --dbfile')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['lookupfn'], mzidtsv_options))

    def get_psms(self):
        """Creates iterator to write to new tsv. Contains input tsv
        lines plus quant data for these."""
        self.header = actions.create_header(self.oldheader)
        self.psms = actions.add_genes_to_psm_table(self.fn, self.oldheader,
                                                   self.lookup)
