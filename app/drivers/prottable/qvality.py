from app.drivers.pycolator.qvality import QvalityDriver
from app.actions.prottable import qvality as preparation
from app.readers import tsv


class ProttableQvalityDriver(QvalityDriver):
    """Runs qvality on two protein tables"""
    outsuffix = '_protqvality.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.score_get_fun = preparation.prepare_qvality_input
        if '--reverse' not in self.qvalityoptions:
            self.qvalityoptions.extend(['--reverse'])
        if self.featuretype not in ['probability']:
            raise Exception('Featuretype (-f) should be proteinprobability.')
        self.score_get_fun = preparation.prepare_qvality_input

    def set_features(self):
        """Creates scorefiles for qvality's target and decoy distributions"""
        targetheader = tsv.get_tsv_header(self.fn)
        self.target = tsv.generate_tsv_proteins(self.fn, targetheader)
        decoyheader = tsv.get_tsv_header(self.decoy)
        self.decoy = tsv.generate_tsv_proteins(self.decoy, decoyheader)
        super().set_features()
