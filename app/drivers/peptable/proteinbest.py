from app.drivers.prottable.base import PepProttableDriver
from app.readers import tsv as tsvreader
import app.actions.peptable.proteinbest as prep


class BestPeptidePerProtein(PepProttableDriver):
    """Filters peptide table to contain only the best scoring peptide for
    each protein in the set.
    """
    outsuffix = '_peptable.tsv'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.scorecol = kwargs.get('scorecol', None)
        self.minlogscore = kwargs.get('logscore', False)

    def initialize_input(self):
        self.oldheader = tsvreader.get_tsv_header(self.fn)
        self.get_column_header_for_number(['scorecol'])

    def create_header(self):
        self.header = self.oldheader[:]

    def set_feature_generator(self):
        self.features = prep.generate_peptides(self.fn, self.oldheader,
                                               self.scorecol, self.minlogscore)
