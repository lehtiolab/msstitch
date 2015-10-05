from app.drivers.prottable.base import ProttableAddData
from app.readers import tsv as tsvreader
import app.actions.prottable.bestpeptide as prep


class BestPeptidePerProtein(ProttableAddData):
    """Filters peptide table to extract the best scoring peptides per protein,
    then adds this score to the protein table.
    """
    outsuffix = '_bestpep.tsv'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.peptable = kwargs.get('pepfile')
        self.scorecol = kwargs.get('scorecolpattern')
        self.minlogscore = kwargs.get('logscore', False)
        self.headertypes = ['bestpepscore']

    def initialize_input(self):
        super().initialize_input()
        self.pepheader = tsvreader.get_tsv_header(self.peptable)
        self.scorecol = tsvreader.get_cols_in_file(self.scorecol,
                                                   self.oldheader, True)

    def set_feature_generator(self):
        self.features = prep.generate_proteins(self.peptable, self.in_proteins,
                                               self.pepheader, self.scorecol,
                                               self.minlogscore)
