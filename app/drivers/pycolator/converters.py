from app.drivers.pycolator import base
from app import converters
from app.writers import pycolator as writers


class Pout2TSVDriver(base.PycolatorDriver):
    outsuffix = '.tsv'

    def __init__(self, **kwargs):
        super(Pout2TSVDriver, self).__init__(**kwargs)
        self.featuretype = kwargs.get('feattype')
        if self.featuretype is None:
            self.featuretype = ['peptide', 'psm']
        elif self.feattype not in ['peptide', 'psm']:
            raise Exception('Featuretype (-f) should be peptide or psm')
        else:
            self.featuretype = [self.featuretype]

    def set_features(self):
        self.featextractors = {
            'peptide': converters.peptide_to_tsv(self.allpeps, self.ns),
            'psm': converters.psm_to_tsv(self.allpsms, self.ns),
        }

    def write(self):
        for feattype in self.featuretype:
            suffix = '_{0}{1}'.format(feattype, self.outsuffix)
            outfn = self.create_outfilepath(self.fn, suffix)
            writers.write_percolator_tsv(self.featextractors[feattype],
                                         feattype, outfn)
