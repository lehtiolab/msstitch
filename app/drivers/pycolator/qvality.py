import subprocess
from app.drivers.pycolator import base
from app.actions.pycolator import qvality as preparation
from app.writers import pycolator as writers


class QvalityDriver(base.PycolatorDriver):
    """Runs qvality from two Percolator XML files. One containing target
    PSMs or peptides, and the other containing decoys."""
    outsuffix = '_qvalityout.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.target = self.fn  # for clarity
        self.decoy = kwargs.get('decoyfn', None)
        if self.decoy is None:
            raise Exception('Not implemented single file qvality yet')
        self.featuretype = kwargs.get('feattype', 'peptide')
        if self.featuretype not in ['peptide', 'psm', 'probability']:
            raise Exception('Featuretype (-f) should be peptide or psm.')
        self.qvalityoptions = []
        options = kwargs.get('options')
        if options is not None:
            for option in kwargs.get('options'):
                self.qvalityoptions.append(option.replace('***', '--'))
        self.score_get_fun = preparation.prepare_qvality_input

    def extract_from_single_target_decoy_file(self):
        """Future method to split into target and decoy file"""
        # FIXME write method, use Splittd driver
        pass

    def set_features(self):
        """Creates scorefiles for qvality's target and decoy distributions"""
        self.scores = {}
        for t_or_d, fn in zip(['target', 'decoy'], [self.target, self.decoy]):
            self.scores[t_or_d] = {}
            self.scores[t_or_d]['scores'] = self.score_get_fun(
                fn, self.featuretype, self.prepare_percolator_output)
            self.scores[t_or_d]['fn'] = self.create_outfilepath(
                t_or_d, self.outsuffix)
            writers.write_qvality_input(self.scores[t_or_d]['scores'],
                                        self.scores[t_or_d]['fn'])

    def write(self):
        """This actually runs the qvality program from PATH."""
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        command = ['qvality']
        command.extend(self.qvalityoptions)
        command.extend([self.scores['target']['fn'], self.scores['decoy']['fn']])
        command.extend(['-o', outfn])
        subprocess.call(command)
