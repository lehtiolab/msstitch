from app.drivers.prottable.base import PepProttableDriver
from app.actions.headers import peptable as head
from app.readers import tsv as tsvreader
import app.actions.peptable.psmtopeptable as prep


class ModelQValuesDriver(PepProttableDriver):
    """Given a peptide table, this uses linear regression to model the
    peptide q-values against a score, e.g. svm-score.
    # FIXME
    It currently also removes the column with PEPs, since it will no
    longer be correct.
    """
    outsuffix = '_qmodel.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.scorecol = kwargs.get('scorecolpattern')
        self.qvalcol = kwargs.get('qcolpattern')

    def initialize_input(self):
        self.oldheader = tsvreader.get_tsv_header(self.fn)
        self.scorecol = tsvreader.get_cols_in_file(self.scorecol,
                                                   self.oldheader, True)
        self.qvalcol = tsvreader.get_cols_in_file(self.qvalcol,
                                                  self.qvalcol, True)

    def create_header(self):
        self.header = head.get_linear_model_header(self.oldheader)

    def set_feature_generator(self):
        self.features = prep.recalculate_qvals_linear_model(self.fn,
                                                            self.oldheader,
                                                            self.scorecol,
                                                            self.qvalcol,
                                                            10e-4)
