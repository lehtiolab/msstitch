import subprocess

from app.drivers.basedrivers import PycolatorDriver
from app.readers import pycolator as readers
from app.writers import pycolator as writers
from app.preparation import pycolator as preparation


class ReassignmentDriver(PycolatorDriver):
    """Reassigns statistics from qvality output on a percolator output file"""
    outsuffix = '_reassigned.xml'

    def __init__(self, **kwargs):
        super(ReassignmentDriver, self).__init__(**kwargs)
        self.qvalityout = kwargs['qvalityout']

    def set_features(self):
        stats = qvality.parse_qvality_output(self.qvalityout)
        self.features = {'peptide': qvality.reassign_elements(self.allpeps,
                                                              stats,
                                                              self.ns),
                         'psm': []
                         }


class QvalityDriver(PycolatorDriver):
    """Runs qvality from two Percolator XML files. One containing target
    PSMs or peptides, and the other containing decoys."""
    outsuffix = '_qvalityout.txt'

    def __init__(self, **kwargs):
        super(QvalityDriver, self).__init__(**kwargs)
        self.target = self.fn  # for clarity
        self.decoy = kwargs.get('decoyfn', None)
        if self.decoy is None:
            raise Exception('Not implemented single file qvality yet')
        self.featuretype = kwargs.get('feattype', 'peptide')
        if self.featuretype not in ['peptide', 'psm']:
            raise Exception('Featuretype (-f) should be peptide or psm.')
        self.qvalityoptions = []
        options = kwargs.get('options')
        if options is not None:
            for option in kwargs.get('options'):
                self.qvalityoptions.append(option.replace('***', '--'))

    def extract_from_single_target_decoy_file(self):
        """Future method to split into target and decoy file"""
        # FIXME write method, use Splittd driver
        pass

    def prepare(self):
        """Qvality running does not need allpsms or static xml, we do all
        of that in set_features instead.
        """
        pass

    def set_features(self):
        """Creates scorefiles for qvality's target and decoy distributions"""
        features_for_qvality = {'target': {}, 'decoy': {}}
        featextractors = {'peptide': readers.generate_peptides,
                          'psm': readers.generate_psms
                          }
        self.scorefiles = {}
        for t_or_d, fn in zip(['target', 'decoy'], [self.target, self.decoy]):
            ns, static_xml = self.prepare_percolator_output(fn)
            features_for_qvality[t_or_d] = featextractors[self.feattype](fn,
                                                                         ns)
            self.scorefiles[t_or_d] = self.create_outfilepath(t_or_d,
                                                              self.outsuffix)
            scores = preparation.get_score(features_for_qvality[t_or_d], ns)
            writers.write_qvality_input(scores, self.scorefiles[t_or_d])

    def write(self):
        """This actually runs the qvality program from PATH."""
        outfn = self.create_outfilepath(self.fn, '_qvalityout.txt')
        command = ['qvality']
        command.extend(self.qvalityoptions)
        command.extend([self.scorefiles['target'], self.scorefiles['decoy']])
        command.extend(['-o', outfn])
        subprocess.call(command)
