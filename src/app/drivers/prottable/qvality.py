from app.drivers.pycolator.qvality import QvalityDriver
from app.actions.prottable import qvality as preparation
from app.actions.prottable import picktdprotein as pickprotein
from app.readers import tsv
from app.drivers.options import prottable_options


class ProttableQvalityDriver(QvalityDriver):
    """Runs qvality on two TSV tables"""
    outsuffix = '_protqvality.txt'
    command = 'qvality'
    commandhelp = ('Run qvality on protein (or tsv) tables '
                   'containing target (-i) proteins and decoy (--decoy) '
                   'proteins. Use with --feattype to use either protein '
                   'error probability (Nesvizhskii 2003) or Q score from '
                   'Savitski 2014 MCP.')

    def set_options(self):
        super().set_options()
        options = self.define_options(['featuretype'], prottable_options)
        self.options.update(options)

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.score_get_fun = preparation.prepare_qvality_input

    def prepare(self):
        """No percolator XML for protein tables"""
        self.targetheader = tsv.get_tsv_header(self.fn)
        self.decoyheader = tsv.get_tsv_header(self.decoyfn)

    def set_features(self):
        """Creates scorefiles for qvality's target and decoy distributions"""
        self.target = tsv.generate_tsv_proteins(self.fn, self.targetheader)
        self.decoy = tsv.generate_tsv_proteins(self.decoyfn, self.decoyheader)
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
    command = 'pickqvality'
    commandhelp = ('Run qvality on protein tables containing target (-i) '
                   'proteins and decoy (--decoy) proteins. Targets and '
                   'decoys will be filtered according to Savitski et al. '
                   '2014 MCP, where the best scoring of a pair of '
                   'target/decoy proteins is retained and the other '
                   'protein discarded. Matching (reversed, tryptic reversed,'
                   ' scrambled) target and decoy FASTA files are needed to '
                   'determine the pairs, use --targetfasta, --decoyfasta. '
                   'Specify gene type with --picktype to distinguish between '
                   'creating target-decoy pairs from '
                   'matched FASTA files or from the protein table input.')

    def set_options(self):
        super().set_options()
        options = self.define_options(['t_fasta', 'd_fasta', 'picktype',
                                       'genefield', 'fastadelim'],
                                      prottable_options)
        self.options.update(options)
        tmp_options = {}
        for clarg, option in self.options.items():
            if option['driverattr'] != 'featuretype':
                tmp_options.update({clarg: option})
        self.options = tmp_options

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.featuretype = 'qvalue'

    def set_features(self):
        """Using this to write picked score tables for qvality"""
        fastadelim, genefield = self.get_fastadelim_genefield(self.fastadelim,
                                                              self.genefield)
        target, decoy = pickprotein.write_pick_td_tables(
            self.target, self.decoy, self.targetheader, self.decoyheader,
            self.t_fasta, self.d_fasta, self.picktype, fastadelim,
            genefield)
        self.scores = {'target': {'fn': target}, 'decoy': {'fn': decoy}}
