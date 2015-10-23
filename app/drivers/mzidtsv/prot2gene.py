from app.actions.mzidtsv import prot2gene as actions
from app.drivers.mzidtsv import MzidTSVDriver


class TSVGeneFromProteinDriver(MzidTSVDriver):
    outsuffix = '_genes.txt'
    lookuptype = 'psm'

    def get_psms(self):
        """Creates iterator to write to new tsv. Contains input tsv
        lines plus quant data for these."""
        self.header = actions.create_header(self.oldheader)
        self.psms = actions.add_genes_to_psm_table(self.fn, self.oldheader,
                                                   self.lookup)
