from app.drivers.prottable.base import PepProttableDriver
from app.actions.headers import peptable as head
from app.readers import tsv as tsvreader
import app.actions.peptable.psmtopeptable as prep


class MzidTSVPeptableDriver(PepProttableDriver):
    """Creates unique peptide table from MzidTSV table. Will not change
    q-values/PEP after filtration, so these should have been calculated
    for the peptides beforehand or rescored afterwards.
    """
    outsuffix = '_peptable.tsv'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.fncol = kwargs.get('speccol', False)
        self.scorecol = kwargs.get('scorecol')
        self.quantcolpattern = kwargs.get('quantcolpattern', None)
        self.precursorquantcolpattern = kwargs.get('precursorquantcolpattern',
                                                   None)

    def initialize_input(self):
        self.oldheader = tsvreader.get_tsv_header(self.fn)
        self.get_column_header_for_number(['fncol', 'scorecol'])
        self.isobfieldmap = prep.get_quantcols(self.quantcolpattern,
                                               self.oldheader, 'isob')
        self.precurquantcol = prep.get_quantcols(self.precursorquantcolpattern,
                                                 self.oldheader, 'precur')

    def create_header(self):
        self.header = head.get_pepquant_header(self.oldheader,
                                               self.isobfieldmap,
                                               self.precurquantcol)

    def set_feature_generator(self):
        self.features = prep.generate_peptides(self.fn, self.oldheader,
                                               self.scorecol,
                                               self.isobfieldmap,
                                               self.precurquantcol, self.fncol)
