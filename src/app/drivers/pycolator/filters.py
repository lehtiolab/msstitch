from app.drivers.pycolator import base
from app.actions.pycolator import filters as preparation
from app.drivers.options import pycolator_options


class FilterPeptideLength(base.PycolatorDriver):
    """Filters on peptide length, to be specified in calling. Outputs to
    multiple files if multiple file input is given. No PSMs will be
    outputted."""
    outsuffix = '_filtlen.xml'
    command = 'filterlen'
    commandhelp = 'Filters out peptides that exceed --maxlen and --minlen'

    def set_options(self):
        super().set_options()
        options = self.define_options(['maxlength', 'minlength'],
                                      pycolator_options)
        self.options.update(options)

    def set_features(self):
        # FIXME psm filter len too!
        self.features = {
            'psm': preparation.filter_peptide_length(
                self.allpsms, 'psm', self.ns, self.minlength, self.maxlength),
            'peptide': preparation.filter_peptide_length(
                self.allpeps, 'pep', self.ns, self.minlength, self.maxlength)
        }


class FilterUniquePeptides(base.PycolatorDriver):
    """This class processes multiple percolator runs from fractions and
    filters out the best scoring peptides."""
    outsuffix = '_filtuni.xml'
    command = 'filteruni'
    commandhelp = ('Only includes best scoring unique peptides in a (merged) '
                   'file')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['score'], pycolator_options))

    def get_all_psms(self):
        """Override parent method so it returns strings instead"""
        return self.get_all_psms_strings()

    def set_features(self):
        uniquepeps = preparation.filter_unique_peptides(self.allpeps,
                                                        self.score,
                                                        self.ns)
        self.features = {'psm': self.allpsms, 'peptide': uniquepeps}


class FilterWholeProteinSequence(base.PycolatorDriver):
    """This class processes multiple percolator runs from fractions and
    filters out first peptides that are found in a specified searchspace. Then
    it keeps the remaining best scoring unique peptides."""
    outsuffix = '_filtprot.xml'
    command = 'filterprot'
    lookuptype = 'searchspace'
    commandhelp = ('Filters out peptides that match to a protein in a FASTA '
                   'file. Needs a lookup to match peptides to quickly.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['deamidate', 'fasta',
                                                 'minlength', 'lookupfn',
                                                 'forcetryp'],
                                                pycolator_options))

    def set_features(self):
        self.features = {
            'peptide': preparation.filter_whole_proteins(self.allpeps,
                                                         self.fasta,
                                                         self.lookup,
                                                         'pep',
                                                         self.ns,
                                                         self.deamidate,
                                                         self.minlength,
                                                         self.forcetryp),
            'psm': preparation.filter_whole_proteins(self.allpsms,
                                                     self.fasta,
                                                     self.lookup,
                                                     'psm',
                                                     self.ns,
                                                     self.deamidate,
                                                     self.minlength,
                                                     self.forcetryp)
        }


class FilterPeptideSequence(base.PycolatorDriver):
    """This class processes multiple percolator runs from fractions and
    filters out first peptides that are found in a specified searchspace. Then
    it keeps the remaining best scoring unique peptides."""
    outsuffix = '_filtseq.xml'
    lookuptype = 'searchspace'
    command = 'filterseq'
    commandhelp = ('Filters out peptides that are found in a certain lookup '
                   'DB')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['falloff', 'deamidate',
                                                 'lookupfn'],
                                                pycolator_options))

    def set_features(self):
        self.features = {
            'peptide': preparation.filter_known_searchspace(self.allpeps,
                                                            'pep',
                                                            self.lookup,
                                                            self.ns,
                                                            self.falloff,
                                                            self.deamidate),
            'psm': preparation.filter_known_searchspace(self.allpsms,
                                                        'psm',
                                                        self.lookup,
                                                        self.ns,
                                                        self.falloff,
                                                        self.deamidate),
        }
