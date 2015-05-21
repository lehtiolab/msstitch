from app.drivers.mzidtsv import MzidTSVDriver
import app.actions.mzidtsv.peptable as prep


class MzidTSVPeptableDriver(MzidTSVDriver):
    """Creates unique peptide table from MzidTSV table. Will not change
    q-values/PEP after filtration, so these should have been calculated
    for the peptides beforehand or rescored afterwards.
    """
    outsuffix = '_peptable.tsv'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.fncol = kwargs.get('fncol', None)
        self.scorecol = kwargs.get('scorecol', None)

    def get_psms(self):
        self.header = prep.get_peptable_header(self.oldheader)
        self.get_column_header_for_number(['fncol', 'scorecol'])
        self.psms = prep.generate_peptides(self.fn, self.oldheader, self.scorecol,
                                           self.fncol)
