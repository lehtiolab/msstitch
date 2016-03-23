from app.drivers.prottable.base import PepProttableDriver
from app.actions.headers import peptable as head
from app.readers import tsv as tsvreader
from app.actions.peptable import psmtopeptable as prep
from app.drivers.options import peptable_options


class MzidTSVPeptableDriver(PepProttableDriver):
    """Creates unique peptide table from MzidTSV table. Will not change
    q-values/PEP after filtration, so these should have been calculated
    for the peptides beforehand or rescored afterwards.
    """
    outsuffix = '_peptable.tsv'
    command = 'psm2pep'
    commandhelp = ('Create peptide table from PSM TSV input, uses best '
                   'scoring PSM for each peptide and strips PSM '
                   'isobaric quant information. Retains MS1 quant info if '
                   'specified, by taking the highest precursor quant value '
                   'for a peptide.')

    def __init__(self):
        super().__init__()
        self.infiletype = 'TSV PSM table (MSGF+)'

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['spectracol',
                                                 'scorecolpattern',
                                                 'quantcolpattern',
                                                 'precursorquantcolpattern'],
                                                peptable_options))
        self.options['--isobquantcolpattern']['help'] += (
            'Isobaric quantification values will be removed from the table '
            'since they represent PSMs')

    def initialize_input(self):
        self.oldheader = tsvreader.get_tsv_header(self.fn)
        self.get_column_header_for_number(['spectracol'])
        self.scorecol = tsvreader.get_cols_in_file(self.scorecolpattern,
                                                   self.oldheader, True)
        self.precurquantcol = prep.get_quantcols(self.precursorquantcolpattern,
                                                 self.oldheader, 'precur')

    def create_header(self):
        self.header = head.get_psm2pep_header(self.oldheader,
                                              self.quantcolpattern,
                                              self.precurquantcol)

    def set_feature_generator(self):
        self.features = prep.generate_peptides(self.fn, self.oldheader,
                                               self.scorecol,
                                               self.precurquantcol,
                                               self.spectracol)
