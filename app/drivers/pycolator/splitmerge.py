from app.drivers.pycolator import base
from app.actions.pycolator import splitmerge as preparation
from app.readers import pycolator as readers


class SplitDriver(base.PycolatorDriver):
    command = 'splittd'
    commandhelp = ('Splits target and decoy data, producing 2 output files')

    def __init__(self):
        super().__init__()
        self.targetsuffix = '_target.xml'
        self.decoysuffix = '_decoy.xml'

    def run(self):
        for filter_type, suffix in zip(['target', 'decoy'],
                                       [self.targetsuffix, self.decoysuffix]):
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
        self.features = preparation.split_target_decoy(elements_to_split,
                                                       self.ns, filter_type)


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

    def set_options(self, options=None):
        super().set_options(['multifiles'])

    def prepare(self):
        self.ns, self.static_xml = self.prepare_percolator_output(self.fn)

    def set_features(self):
        """"Merge all psms and peptides"""
        allpsms_str = readers.generate_psms_multiple_fractions_strings(
            self.mergefiles, self.ns)
        allpeps = preparation.merge_peptides(self.mergefiles, self.ns)
        self.features = {'psm': allpsms_str, 'peptide': allpeps}
