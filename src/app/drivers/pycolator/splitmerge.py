from app.drivers.pycolator import base
from app.actions.pycolator import splitmerge as preparation
from app.readers import pycolator as readers
from app.drivers.options import pycolator_options


class SplitDriver(base.PycolatorDriver):
    command = 'splittd'
    commandhelp = ('Splits target and decoy data, producing 2 output files')

    def run(self):
        for filter_type, suffix in self.filter_types:
            self.prepare()
            self.set_features(filter_type)
            self.outsuffix = suffix
            self.write()
        self.finish()

    def set_features(self, filter_type):
        """Calls splitter to split percolator output into target/decoy
        elements.
        Writes two new xml files with features. Currently only psms and
        peptides. Proteins not here, since one cannot do protein inference
        before having merged and remapped multifraction data anyway.
        """
        elements_to_split = {'psm': self.allpsms, 'peptide': self.allpeps}
        self.features = self.splitfunc(elements_to_split, self.ns, filter_type)


class SplitTDDriver(SplitDriver):
    def __init__(self):
        super().__init__()
        self.targetsuffix = '_target.xml'
        self.decoysuffix = '_decoy.xml'
        self.splitfunc = preparation.split_target_decoy
        self.filter_types = [('target', self.targetsuffix),
                             ('decoy', self.decoysuffix)]


class MergeDriver(base.PycolatorDriver):
    """Base class for merging multiple percolator fractions under different
    sorts of filtering. It writes a single percolator out xml from
    multiple fractions.
    Namespace and static xml come from first percolator file.
    Make sure fractions are from same percolator run."""
    outsuffix = '_merged.xml'
    command = 'merge'
    commandhelp = 'Merges percolator xml files, nothing else.'

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.mergefiles = self.fn[:]
        self.fn = self.fn[0]

    def set_options(self):
        super().set_options()
        options = self.define_options(['multifiles'], pycolator_options)
        self.options.update(options)

    def prepare(self):
        self.ns, self.static_xml = self.prepare_percolator_output(self.fn)

    def set_features(self):
        """"Merge all psms and peptides"""
        allpsms_str = readers.generate_psms_multiple_fractions_strings(
            self.mergefiles, self.ns)
        allpeps = preparation.merge_peptides(self.mergefiles, self.ns)
        self.features = {'psm': allpsms_str, 'peptide': allpeps}
