from app.drivers.pycolator import base
from app.actions.pycolator import reassign


class ReassignmentDriver(base.PycolatorDriver):
    """Reassigns statistics from qvality output on a percolator output file"""
    outsuffix = '_reassigned.xml'
    command = 'reassign'
    commandhelp = ('Reassigns statistics from a qvality output file onto a '
                   'single percolator input file.')

    def set_options(self, options=None):
        super().set_options(['qvalityout', 'featuretype'])

    def set_features(self):
        stats = reassign.parse_qvality_output(self.qvalityout)
        if self.featuretype == 'psm':
            self.features = {
                'psm': reassign.reassign_elements(self.allpsms, stats,
                                                  self.ns),
                'peptide': self.get_all_peptides_strings(),
            }
        elif self.featuretype == 'peptide':
            self.features = {
                'peptide': reassign.reassign_elements(self.allpeps, stats,
                                                      self.ns),
                'psm': self.get_all_psms_strings(),
            }
