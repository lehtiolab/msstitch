from app.drivers.prottable.base import PepProttableDriver
from app.actions.headers import peptable as head
from app.readers import tsv as tsvreader
from app.actions.peptable import model_qvals as prep
from app.drivers.options import peptable_options


class ModelQValuesDriver(PepProttableDriver):
    """Given a peptide table, this uses linear regression to model the
    peptide q-values against a score, e.g. svm-score.
    # FIXME
    It currently also removes the column with PEPs, since it will no
    longer be correct.
    """
    outsuffix = '_qmodel.txt'
    command = 'modelqvals'
    commandhelp = ('Recalculate peptide q-values by creating a linear model '
                   'of them against a score (partial least squares '
                   'regression).')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['scorecolpattern',
                                                 'fdrcolpattern'],
                                                peptable_options))
        self.options['--fdrcolpattern'].update({'required': True})

    def initialize_input(self):
        self.oldheader = tsvreader.get_tsv_header(self.fn)
        self.scorecol = tsvreader.get_cols_in_file(self.scorecolpattern,
                                                   self.oldheader, True)
        self.qvalcol = tsvreader.get_cols_in_file(self.fdrcolpattern,
                                                  self.oldheader, True)

    def create_header(self):
        self.header = head.get_linear_model_header(self.oldheader)

    def set_feature_generator(self):
        self.features = prep.recalculate_qvals_linear_model(self.fn,
                                                            self.scorecol,
                                                            self.qvalcol,
                                                            10e-4)
