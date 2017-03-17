import subprocess
from app.drivers.pycolator import base
from app.actions.pycolator import qvality as preparation
from app.writers import pycolator as writers
from app.drivers.options import pycolator_options


class QvalityDriver(base.PycolatorDriver):
    """Runs qvality from two Percolator XML files. One containing target
    PSMs or peptides, and the other containing decoys."""
    outsuffix = '_qvalityout.txt'
    command = 'qvality'
    commandhelp = ('Runs qvality on an inputfile: target and decoy data. '
                   'When using separate files for target and decoy, '
                   'use --decoy to specify the decoy input file, and -f '
                   'to specify feature type (psm or peptide).')

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.score_get_fun = preparation.prepare_qvality_input
        self.target = self.fn  # for clarity
        if self.decoyfn is None:
            raise Exception('Not implemented single file qvality yet')
        self.decoy = self.decoyfn
        self.qvalityoptions = []
        if self.qoptions is not None:
            for option in self.qoptions:
                self.qvalityoptions.append(option.replace('***', '--'))

    def set_options(self):
        super().set_options()
        options = self.define_options(['decoyfn', 'featuretype', 'qoptions'],
                                      pycolator_options)
        self.options.update(options)

    def extract_from_single_target_decoy_file(self):
        """Future method to split into target and decoy file"""
        # FIXME write method, use Splittd driver
        pass

    def set_features(self):
        """Creates scorefiles for qvality's target and decoy distributions"""
        self.scores = {}
        for t_or_d, feats in zip(['target', 'decoy'], [self.target,
                                                       self.decoy]):
            self.scores[t_or_d] = {}
            self.scores[t_or_d]['scores'] = self.score_get_fun(
                feats, self.featuretype, self.prepare_percolator_output)
            self.scores[t_or_d]['fn'] = '{}_qvality_input.txt'.format(t_or_d)
            writers.write_qvality_input(self.scores[t_or_d]['scores'],
                                        self.scores[t_or_d]['fn'])

    def write(self):
        """This actually runs the qvality program from PATH."""
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        command = ['qvality']
        command.extend(self.qvalityoptions)
        command.extend([self.scores['target']['fn'], self.scores['decoy']['fn'],
                        '-o', outfn])
        subprocess.call(command)
