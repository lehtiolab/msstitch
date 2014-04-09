from app.drivers.basedrivers import PycolatorDriver
from app.preparation import pycolator as preparation


class FilterPeptideLength(PycolatorDriver):
    """Filters on peptide length, to be specified in calling. Outputs to
    multiple files if multiple file input is given. No PSMs will be
    outputted."""
    outsuffix = '_filt_len.xml'

    def __init__(self, **kwargs):
        super(FilterPeptideLength, self).__init__(**kwargs)
        self.minlength = kwargs.get('minlength', 0)
        self.maxlength = kwargs.get('maxlength', None)

    def set_features(self):
        # FIXME psm filter len too!
        self.features = {'psm': [],
                         'peptide': preparation.filter_peptide_length(
                             self.allpeps, self.ns,
                             self.minlength, self.maxlength)
                         }


class FilterUniquePeptides(PycolatorDriver):
    """This class processes multiple percolator runs from fractions and
    filters out the best scoring peptides."""
    outsuffix = '_filtuniq.xml'

    def __init__(self, **kwargs):
        super(FilterUniquePeptides, self).__init__(**kwargs)
        self.score = kwargs.get('score')
        if self.score is None:
            self.score = 'svm'

    def set_features(self):
        uniquepeps = preparation.filter_unique_peptides(self.allpeps,
                                                        self.score,
                                                        self.ns)
        # FIXME stringify psms
        self.features = {'psm': self.allpsms_str, 'peptide': uniquepeps}


class FilterKnownPeptides(PycolatorDriver):
    """This class processes multiple percolator runs from fractions and
    filters out first peptides that are found in a specified searchspace. Then
    it keeps the remaining best scoring unique peptides."""
    outsuffix = '_filtknown.xml'

    def __init__(self, **kwargs):
        super(PycolatorDriver, self).__init__(**kwargs)
        self.db = kwargs.get('database')
        self.falloff = kwargs.get('falloff')

    def set_features(self):
        novelpeps = preparation.filter_known_searchspace(self.allpeps,
                                                         self.db,
                                                         self.ns,
                                                         self.falloff)
        self.features = {'psm': [], 'peptide': novelpeps}
