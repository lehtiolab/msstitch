from app.drivers.prottable.base import ProttableAddData
from app.readers import tsv as tsvreader
from app.actions.prottable import bestpeptide as prep
from app.drivers.options import prottable_options


class BestPeptidePerProtein(ProttableAddData):
    """Filters peptide table to extract the best scoring peptides per protein,
    then adds this score to the protein table.
    """
    outsuffix = '_bestpep.tsv'
    command = 'bestpeptide'
    commandhelp = ('Given the protein table and corresponding peptide table '
                   ', fetch the best scoring peptide for each protein '
                   'and annotate that score in the protein table.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['pepfile', 'proteincol',
                                                 'scorecolpattern',
                                                 'pcolpattern', 'minlogscore'],
                                                prottable_options))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.headertypes = ['bestpepscore']

    def initialize_input(self):
        super().initialize_input()
        self.pepheader = tsvreader.get_tsv_header(self.pepfile)
        if self.proteincol:
            self.get_column_header_for_number(['proteincol'], self.pepheader)
        elif self.pcolpattern:
            self.proteincol = tsvreader.get_cols_in_file(self.pcolpattern,
                                                         self.pepheader, True)
        self.scorecol = tsvreader.get_cols_in_file(self.scorecolpattern,
                                                   self.pepheader, True)

    def set_feature_generator(self):
        self.features = prep.generate_proteins(self.pepfile, self.in_proteins,
                                               self.pepheader, self.scorecol,
                                               self.minlogscore,
                                               protcol=self.proteincol)
