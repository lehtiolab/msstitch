from app.drivers.pycolator import base
from app.actions import oldpycolator as preparation
from app.readers import pycolator as readers


class SplitDriver(base.PycolatorDriver):
    def __init__(self, **kwargs):
        super(SplitDriver, self).__init__(**kwargs)
        self.targetsuffix = kwargs.get('targetsuffix', '_target.xml')
        self.decoysuffix = kwargs.get('decoysuffix', '_decoy.xml')

    def run(self):
        td = {'target': self.targetsuffix, 'decoy': self.decoysuffix}
        for filter_type in ['target', 'decoy']:
            self.prepare()
            self.set_features(filter_type)
            self.outsuffix = td[filter_type]
            self.write()
        self.finish()

    def set_features(self, filter_type):
        """ Calls splitter to split percolator output into target/decoy elements.
            Writes two new xml files with features. Currently only psms and
            peptides. Proteins not here, since one cannot do protein inference
            before having merged and remapped multifraction data anyway.
        """
        elements_to_split = {'psm': self.allpsms, 'peptide': self.allpeps}
        self.features = preparation.split_target_decoy(elements_to_split,
                                                       self.ns, filter_type)


class MergeDriver(base.PycolatorDriver):
    """Base class for merging multiple percolator fractions under different
    sorts of filtering. It writes a single percolator out xml from multiple fractions.
    Namespace and static xml come from first percolator file.
    Make sure fractions are from same percolator run."""
    outsuffix = '_merged.xml'

    def __init__(self, **kwargs):
        super(MergeDriver, self).__init__(**kwargs)
        self.mergefiles = [self.fn]
        self.mergefiles.extend(kwargs.get('multifile_input', None))

    def prepare(self):
        self.ns, self.static_xml = self.prepare_percolator_output(self.fn)

    def set_features(self):
        """"Merge all psms and peptides"""
        allpsms_str = readers.generate_psms_multiple_fractions_strings(
            self.mergefiles, self.ns)
        allpeps = preparation.merge_peptides(self.mergefiles, self.ns)
        self.features = {'psm': allpsms_str, 'peptide': allpeps}
