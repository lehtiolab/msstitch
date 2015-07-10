from app.drivers.pycolator import base
from app.actions.pycolator import reassign


class ReassignmentDriver(base.PycolatorDriver):
    """Reassigns statistics from qvality output on a percolator output file"""
    outsuffix = '_reassigned.xml'

    def __init__(self, **kwargs):
        super(ReassignmentDriver, self).__init__(**kwargs)
        self.qvalityout = kwargs['qvalityout']
        self.reassign_feats = kwargs['feattype']

    def set_features(self):
        stats = reassign.parse_qvality_output(self.qvalityout)
        if self.reassign_feats == 'psm':
            self.features = {
                'psm': reassign.reassign_elements(self.allpsms, stats,
                                                  self.ns),
                'peptide': self.get_all_peptides_strings(),
            }
        elif self.reassign_feats == 'peptide':
            self.features = {
                'peptide': reassign.reassign_elements(self.allpeps, stats,
                                                      self.ns),
                'psm': self.get_all_psms_strings(),
            }
