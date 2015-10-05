from app.drivers.pycolator.qvality import QvalityDriver
from app.actions.prottable import qvality as preparation
from app.actions.prottable import picktdprotein as pickprotein
from app.readers import tsv


class ProttableQvalityDriver(QvalityDriver):
    """Runs qvality on two TSV tables"""
    outsuffix = '_protqvality.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.featuretype not in ['probability', 'qvalue']:
            raise Exception('Featuretype (-f) should be (protein) probability '
                            'or Q-score (qvalue).')
        self.score_get_fun = preparation.prepare_qvality_input

    def prepare(self):
        """No percolator XML for protein tables"""
        pass

    def set_features(self):
        """Creates scorefiles for qvality's target and decoy distributions"""
        targetheader = tsv.get_tsv_header(self.fn)
        self.target = tsv.generate_tsv_proteins(self.fn, targetheader)
        decoyheader = tsv.get_tsv_header(self.decoy)
        self.decoy = tsv.generate_tsv_proteins(self.decoy, decoyheader)
        super().set_features()


class PickedQvalityDriver(ProttableQvalityDriver):
    """Given target and decoy protein tables, and matching target and decoy
    FASTA files, this produces a target and decoy protein table with only
    a column for score and protein accession. It picks the best scoring
    protein for each target/decoy pair and outputs that to its corresponding
    new table. Score is currently assumed to be Q score. After the picking,
    qvality is run to output an FDR score table.
    """
    outsuffix = '_pickedqvality.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.featuretype not in ['qvalue']:
            raise Exception('Featuretype (-f) should be '
                            'Q-score (qvalue).')
        self.t_fasta = kwargs.get('targetfasta')
        self.d_fasta = kwargs.get('decoyfasta')

    def prepare(self):
        """Using this to write picked protein tables"""
        self.target, self.decoy = pickprotein.write_pick_td_tables(
            self.target, self.decoy, self.t_fasta, self.d_fasta)
