import os

import readers
import filtering
import writers
import databases

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

    def create_outfilepath(self, fn, suffix=None):
        basefn = os.path.basename(fn)
        outfn = basefn + suffix
        return os.path.join(self.outdir, outfn)


def merge_multiple_fractions(fns):
    """Performs the work to merge parallelized percolator fractions.
    Target/decoy split, filtering unique peptides, running qvality on resulting
    score distributions for psms and peptides and setting values."""
    pass

class SplitDriver(BaseDriver):
    def __init__(self, **kwargs):
        super(MergeDriver, self).__init__(fns, outdir, **kwargs)
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
            namespace, static_xml = self.prepare_percolator_output(fn)
            split_elements = filtering.split_target_decoy(fn, namespace)
            targetfn = self.create_outfilepath(fn, self.outdir, self.targetsuffix)
            decoyfn = self.create_outfilepath(fn, self.outdir, self.decoysuffix)
            writers.write_percolator_xml(static_xml, split_elements['target'], 
                                            targetfn)
            writers.write_percolator_xml(static_xml, split_elements['decoy'],
                                            decoyfn)



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
        self.db = kwargs.get('database', False)

    def run(self):
        self.prepare_merge()
        self.merge()
        self.write()

    def prepare_merge(self):
        self.ns, self.static_xml = self.prepare_percolator_output(self.fns[0])
        self.allpsms = readers.generate_psms_multiple_fractions(self.fns, self.ns)
        self.allpeps = readers.generate_peptides_multiple_fractions(self.fns, self.ns)

    def write(self):
        merged_fn = self.create_outfilepath(self.fns[0], self.outsuffix)
        writers.write_percolator_xml(self.static_xml, self.features, merged_fn)


## FIXME we should split up merging and filtering step in case someone needs
# only filtering, (which we dont have a current use case for).

class MergeUniqueAndFilterKnownPeptides(MergeDriver):
    """This class processes multiple percolator runs from fractions and
    filters out first peptides that are found in a specified searchspace. Then
    it keeps the remaining best scoring unique peptides."""
    def run(self):
        print 'Digesting database into memory to get known search space'
        self.searchspace = databases.get_searchspace(self.db)
        print 'Filtering and merging'
        super(MergeUniqueAndFilterKnownPeptides, self).run()
        assert self.db not in [False, None]

    def merge(self):
        newpeps = filtering.filter_known_searchspace(self.allpeps,
                                            self.searchspace, self.ns)
        uniquepeps = filtering.filter_unique_peptides(newpeps, self.score,
                                                    self.ns)
        self.features = {'psm': self.allpsms, 'peptide': uniquepeps}


class MergeUniquePeptides(MergeDriver):
    """This class processes multiple percolator runs from fractions and
    filters out the best scoring peptides."""
    def merge(self):
        uniquepeps = filtering.filter_unique_peptides(self.allpeps, self.score,
                                                        self.ns)
        self.features = {'psm': self.allpsms, 'peptide': uniquepeps}

