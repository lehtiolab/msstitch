import os
import subprocess

from app.readers import pycolator as readers
from app.preparation import pycolator as preparation
from app.writers import pycolator as writers
from app import qvality
from app.lookups import sequences

class BaseDriver(object):
    def __init__(self, **kwargs):
        self.fns = kwargs['infile']
        self.outdir = kwargs['outdir']
        self.outsuffix = kwargs.get('outsuffix', '.xml')

    def prepare_percolator_output(self, fn):
        """Returns namespace and static xml from percolator output file"""
        ns = readers.get_namespace(fn)
        static = readers.get_percolator_static_xml(fn, ns)
        return ns, static

    def get_all_peptides(self):
        return readers.generate_peptides_multiple_fractions(self.fns, self.ns)

    def get_all_psms(self):
        return readers.generate_psms_multiple_fractions(self.fns, self.ns)

    def create_outfilepath(self, fn, suffix=None):
        basefn = os.path.basename(fn)
        outfn = basefn + suffix
        return os.path.join(self.outdir, outfn)


class QvalityDriver(BaseDriver):
    def __init__(self, **kwargs):
        super(QvalityDriver, self).__init__(**kwargs)
        self.featuretype = kwargs.get('feattype', 'peptide')
        if self.featuretype not in ['peptide', 'psm']:
            raise Exception('Featuretype (-f) should be peptide or psm.')
        self.qvalityoptions = []
        options = kwargs.get('options')
        if options is not None:
            for option in kwargs.get('options'):
                self.qvalityoptions.append(option.replace('***', '--'))

    def run(self):
        scorefiles = self.create_input_scorefiles()
        self.run_qvality(scorefiles)

    def create_input_scorefiles(self):
        scorefiles = {}
        features_for_qvality = {'target': {}, 'decoy': {}}
        for t_or_d, fn in zip(['target', 'decoy'], self.fns):
            ns, static_xml = self.prepare_percolator_output(fn)
            features_for_qvality[t_or_d]['peptide'] = \
                        readers.generate_peptides_multiple_fractions([fn], ns)
            features_for_qvality[t_or_d]['psm'] = \
                        readers.generate_psms_multiple_fractions([fn], ns)
            scorefiles[t_or_d] = self.create_outfilepath(t_or_d, self.outsuffix)
            feats_to_write = preparation.get_score(
                        features_for_qvality[t_or_d][self.featuretype], ns)
            writers.write_qvality_input(feats_to_write, scorefiles[t_or_d])
        return scorefiles

    def run_qvality(self, scorefiles):
        outfn = self.create_outfilepath(self.fns[0], '_qvalityout.txt')
        command = ['qvality']
        command.extend(self.qvalityoptions)
        command.extend([scorefiles['target'], scorefiles['decoy']])
        command.extend(['-o', outfn])
        subprocess.call(command)

class SplitDriver(BaseDriver):
    def __init__(self, **kwargs):
        super(SplitDriver, self).__init__(**kwargs)
        self.targetsuffix = kwargs.get('targetsuffix', '_target.xml')
        self.decoysuffix = kwargs.get('decoysuffix', '_decoy.xml')

    def run(self):
        self.split()

    def split(self): #targetsuffix='_target.xml', decoysuffix='_decoy.xml'):
        """ Calls splitter to split percolator output into target/decoy elements.
            Writes two new xml files with features. Currently only psms and
            peptides. Proteins not here, since one cannot do protein inference
            before having merged and remapped multifraction data anyway.
        """
        for fn in self.fns:
            ns, static_xml = self.prepare_percolator_output(fn)
            pep_target = readers.generate_peptides_multiple_fractions([fn], ns)
            pep_decoy = readers.generate_peptides_multiple_fractions([fn], ns)
            psm_target = readers.generate_psms_multiple_fractions([fn], ns)
            psm_decoy = readers.generate_psms_multiple_fractions([fn], ns)
            elements_to_split = {
            'target':   {'psm': psm_target, 'peptide': pep_target},
            'decoy' :   {'psm': psm_decoy,  'peptide': pep_decoy},
            }
            split_els = preparation.split_target_decoy(elements_to_split, ns)
            targetfn = self.create_outfilepath(fn, self.targetsuffix)
            decoyfn = self.create_outfilepath(fn, self.decoysuffix)
            writers.write_percolator_xml(static_xml, split_els['target'],
                                            targetfn)
            writers.write_percolator_xml(static_xml, split_els['decoy'],
                                            decoyfn)


class ReassignmentDriver(BaseDriver):
    """Reassigns statistics from qvality output on a percolator output file"""
    def __init__(self, **kwargs):
        super(ReassignmentDriver, self).__init__(**kwargs)
        self.qvalityout = kwargs['qvalityout']
        self.outsuffix = kwargs.get('outsuffix', '_reassigned.xml')

    def run(self):
        self.reassign()

    def reassign(self):
        ns, static_xml = self.prepare_percolator_output(self.fns[0])
        allpeps = readers.generate_peptides_multiple_fractions(self.fns, ns)
        stats = qvality.parse_qvality_output(self.qvalityout)
        features = {'peptide': qvality.reassign_elements(allpeps,
                                                         stats,
                                                         ns),
                    'psm': []
                   }
        outfn = self.create_outfilepath(self.fns[0], self.outsuffix)
        writers.write_percolator_xml(static_xml, features, outfn)


class FilterPeptideLengthDriver(BaseDriver):
    """Filters on peptide length, to be specified in calling. Outputs to
    multiple files if multiple file input is given. No PSMs will be
    outputted."""
    def __init__(self, **kwargs):
        super(FilterPeptideLengthDriver, self).__init__(**kwargs)
        self.outsuffix = kwargs.get('outsuffix', 'filtered.xml')
        self.minlength = kwargs.get('minlength', 0)
        self.maxlength = kwargs.get('maxlength', None)

    def prepare(self, fn):
        self.ns, self.static_xml = self.prepare_percolator_output(fn)
        self.allpeps = self.get_all_peptides()
        self.allpsms = self.get_all_psms()

    def run(self):
        for fn in self.fns:
            self.prepare(fn)
            self.features = {'psm': [],
                             'peptide': preparation.filter_peptide_length(
                                 self.allpeps, self.ns,
                                 self.minlength, self.maxlength)
                             }
            self.write(fn)

    def write(self, fn):
        outfn = self.create_outfilepath(fn, self.outsuffix)
        writers.write_percolator_xml(self.static_xml, self.features, outfn)


class MergeDriver(BaseDriver):
    """Base class for merging multiple percolator fractions under different
    sorts of filtering. It writes a single percolator out xml from multiple fractions.
    Namespace and static xml come from first percolator file.
    Make sure fractions are from same percolator run."""

    def __init__(self, **kwargs):
        super(MergeDriver, self).__init__(**kwargs)
        self.outsuffix = kwargs.get('outsuffix', '_merged.xml')
        self.score = kwargs.get('score')
        if self.score is None:
            self.score = 'svm'

    def run(self):
        self.prepare_merge()
        self.merge()
        self.write()

    def prepare_merge(self):
        self.ns, self.static_xml = self.prepare_percolator_output(self.fns[0])
        self.allpeps = self.get_all_peptides()
        self.allpsms = self.get_all_psms()

        self.allpsms_str = readers.generate_psms_multiple_fractions_strings(
            self.fns, self.ns)
        self.allpeps_str = readers.generate_peptides_multiple_fractions_strings(
            self.fns, self.ns)

    def merge(self):
        """"Merge all psms and peptides"""
        self.features = {'psm': self.allpsms_str, 'peptide': self.allpeps_str}

    def write(self):
        merged_fn = self.create_outfilepath(self.fns[0], self.outsuffix)
        writers.write_percolator_xml(self.static_xml, self.features, merged_fn)


## FIXME we should split up merging and filtering step in case someone needs
# only filtering, (which we dont have a current use case for).

class MergeUniqueAndFilterKnownPeptides(MergeDriver):
    """This class processes multiple percolator runs from fractions and
    filters out first peptides that are found in a specified searchspace. Then
    it keeps the remaining best scoring unique peptides."""
    # FIXME create merge drivers, filter drivers, split drivers in diff files
    # app/drivers/pycolator/split, merge, filter, etc
    def __init__(self, **kwargs):
        super(MergeDriver, self).__init__(**kwargs)
        self.db = kwargs.get('database', False)
        self.outsuffix = kwargs.get('outsuffix', '_filterknown.xml')
        self.prolinecut = kwargs.get('proline')
        self.falloff = kwargs.get('falloff')

    def run(self):
        self.searchspace = sequences.create_searchspace(self.db,
                                                        self.prolinecut,
                                                        self.falloff)
        super(MergeUniqueAndFilterKnownPeptides, self).run()
        assert self.db not in [False, None]

    def merge(self):
        novelpeps = preparation.filter_known_searchspace(self.allpeps,
                                                         self.searchspace,
                                                         self.ns,
                                                         self.falloff)
        self.features = {'psm': [], 'peptide': novelpeps}


class MergeUniquePeptides(MergeDriver):
    """This class processes multiple percolator runs from fractions and
    filters out the best scoring peptides."""
    def merge(self):
        uniquepeps = preparation.filter_unique_peptides(self.allpeps, self.score,
                                                        self.ns)
        self.features = {'psm': self.allpsms_str, 'peptide': uniquepeps}

