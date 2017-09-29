from app.drivers.prottable.base import ProttableAddData
from app.actions.prottable import fdr as action
from app.readers import tsv as reader
from app.drivers.options import prottable_options


class ProttableFDRDriver(ProttableAddData):
    """Runs qvality on two TSV tables"""
    outsuffix = '_protfdr.txt'
    command = 'protfdr'
    commandhelp = ('Calculate FDR of protein (or tsv) tables '
                   'containing target (-i) proteins and decoy (--decoy) '
                   'proteins. Uses Q-score from Savitski 2014 MCP as a '
                   'ranking mechanism, q-values are #D/(#D+#T) without '
                   'correction')

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.headertypes = ['proteinfdr']

    def initialize_input(self):
        self.target = self.fn
        self.oldheader = reader.get_tsv_header(self.target)
        self.targetheader = reader.get_tsv_header(self.target)
        self.decoyheader = reader.get_tsv_header(self.decoyfn)

    def set_options(self):
        super().set_options()
        options = self.define_options(['decoyfn'], prottable_options)
        self.options.update(options)

    def prepare(self):
        """No percolator XML for protein tables"""
        self.target = self.fn
        self.targetheader = reader.get_tsv_header(self.target)
        self.decoyheader = reader.get_tsv_header(self.decoyfn)

    def set_feature_generator(self):
        self.features = action.generate_protein_fdr(
            self.target, self.decoyfn, self.targetheader, self.decoyheader,
            self.headerfields)


class PickedFDRDriver(ProttableFDRDriver):
    outsuffix = '_pickedfdr.txt'
    command = 'pickedfdr'
    commandhelp = ('Calculate FDR on gene tables containing target (-i) '
                   'proteins and decoy (--decoy) genes. Targets and '
                   'decoys will be filtered according to Savitski et al. '
                   '2014 MCP, where the best scoring of a pair of '
                   'target/decoy genes is retained and the other '
                   'gene discarded. Matching (reversed, tryptic reversed,'
                   ' scrambled) target and decoy FASTA files are needed to '
                   'determine the pairs, use --targetfasta, --decoyfasta. '
                   'Specify gene type with --picktype to distinguish between '
                   'creating target-decoy pairs from '
                   'matched FASTA files or from the gene table input. '
                   'Uses Q-score as a ranking mechanism, q-values are '
                   '#D/(#D+#T) without correction')

    def set_options(self):
        super().set_options()
        options = self.define_options(['t_fasta', 'd_fasta', 'picktype',
                                       'genefield', 'fastadelim'],
                                      prottable_options)
        self.options.update(options)

    def set_feature_generator(self):
        fastadelim, genefield = self.get_fastadelim_genefield(self.fastadelim,
                                                              self.genefield)
        self.features = action.generate_pick_fdr(
            self.target, self.decoyfn, self.targetheader, self.decoyheader,
            self.headerfields, self.t_fasta, self.d_fasta, self.picktype,
            fastadelim, genefield)


class ProttableAddFDRDriver(ProttableAddData):
    """Assigns FDR to protein table based on qvality output and protein table
    probabilities"""
    outsuffix = '_protfdr.txt'
    command = 'fdr'
    commandhelp = ('Add protein FDR to protein table by comparing '
                   'score (peptide q-value, protein probability, etc) '
                   'with qvality lookup table. Needs --scorecolpattern and '
                   'qvality output file specified with --qvality')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['qvalityout',
                                                 'scorecolpattern'],
                                                prottable_options))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.headertypes = ['proteinfdr', 'proteinpep']

    def initialize_input(self):
        super().initialize_input()
        self.scorecol = reader.get_cols_in_file(self.scorecolpattern,
                                                self.oldheader, True)

    def set_feature_generator(self):
        self.features = action.add_protein_fdr(self.qvalityout,
                                               self.in_proteins,
                                               self.headerfields,
                                               self.scorecol)
