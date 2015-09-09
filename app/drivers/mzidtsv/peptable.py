from app.drivers.prottable import PepProttableDriver
from app.actions.headers import peptable as head
import app.actions.mzidtsv.peptable as prep


class MzidTSVPeptableDriver(PepProttableDriver):
    """Creates unique peptide table from MzidTSV table. Will not change
    q-values/PEP after filtration, so these should have been calculated
    for the peptides beforehand or rescored afterwards.
    """
    outsuffix = '_peptable.tsv'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.fncol = kwargs.get('fncol', None)
        self.scorecol = kwargs.get('scorecol', None)
        self.quantcolpattern = kwargs.get('quantcolpattern', None)
        self.precursorquantcolpattern = kwargs.get('precursorquantcolpattern',
                                                   None)

    def initialize_input(self):
        self.get_column_header_for_number(['fncol', 'scorecol'])
        self.isobfieldmap = prep.get_quantcols(self.quantcolpattern,
                                               self.oldheader, 'isob')
        self.precurquantcol = prep.get_quantcols(self.precursorquantcolpattern,
                                                 self.oldheader, 'precur')

    def create_header(self):
        self.header = head.get_pepquant_header(self.oldheader,
                                               self.isobfieldmap,
                                               self.precurquantcol)

    def get_column_header_for_number(self, column_var_names):
        """This function subtracts 1 from inputted column number to comply
        with programmers counting (i.e. from 0, not from 1). Could possibly
        made more TSV general"""
        convert_columns = ['confcol', 'spec_column']
        convert_columns.extend(column_var_names)
        for col in convert_columns:
            value = getattr(self, col)
            if value is None:
                continue
            setattr(self, col, self.oldheader[int(value) - 1])

    def get_psms(self):
        self.psms = prep.generate_peptides(self.fn, self.oldheader,
                                           self.scorecol, self.isobfieldmap,
                                           self.precurquantcol, self.fncol)
