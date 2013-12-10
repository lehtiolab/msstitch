import os

import readers
import filtering
import writers


class BaseDriver(object):
    def __init__(self, fns, outdir, **kwargs):
        self.fns = fns
        self.outdir = outdir
        if 'outsuffix' in kwargs:
            self.outsuffix = kwargs['outsuffix']
        if 'targetsuffix' in kwargs:
            self.targetsuffix = kwargs['targetsuffix']
        if 'decoysuffix' in kwargs:
            self.decoysuffix = kwargs['decoysuffix']

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
    def split_target_decoy(self, fns, outdir, targetsuffix='_target.xml', decoysuffix='_decoy.xml'):
        """ Calls splitter to split percolator output into target/decoy elements.
            Writes two new xml files with features. Currently only psms and
            peptides. Proteins not here, since one cannot do protein inference
            before having merged and remapped multifraction data anyway.
        """
        for fn in fns:
            namespace, static_xml = self.prepare_percolator_output(fn)
            split_elements = filtering.split_target_decoy(fn, namespace)
            targetfn = create_outfilepath(fn, outdir, targetsuffix)
            decoyfn = create_outfilepath(fn, outdir, decoysuffix)
            writers.write_percolator_xml(static_xml, split_elements['target'], 
                                            targetfn)
            writers.write_percolator_xml(static_xml, split_elements['decoy'],
                                            decoyfn)



class MergeDriver(BaseDriver):
    """Base class for merging multiple percolator fractions under different
    sorts of filtering"""
    def __init__(self, fns, outdir, **kwargs):
        super(MergeDriver, self).__init__(fns, outdir, kwargs)
        if 'score' in kwargs:
            self.score = kwargs['score']
        else:
            self.score = 'svm'

        if 'searchspace' in kwargs:
            self.searchspace = kwargs['searchspace']

    def run(self):
        self.prepare_merge()
        self.merge()
        self.write()

    def prepare_merge(self):
        self.ns, self.static_xml = self.prepare_percolator_output(self.fns[0])
        self.allpsms = readers.generate_psms_multiple_fractions(self.fns, self.ns)
        self.allpeps = readers.generate_psms_multiple_fractions(self.fns, self.ns)

    def write(self):
        merged_fn = self.create_outfilepath(self.fns[0], self.outsuffix)
        writers.write_percolator_xml(self.static_xml, self.features, merged_fn)


class MergeUniqueAndFilterKnownPeptides(MergeDriver):
    def merge(self):
        self.prepare_merge()
        newpeps = filtering.filter_known_searchspace(self.allpeps,
                                            self.searchspace, self.ns)
        uniquepeps = filtering.filter_unique_peptides(newpeps, self.score,
                                                    self.ns)
        self.features = {'psm': self.allpsms, 'peptide': uniquepeps}


class MergeUniquePeptides(MergeDriver):
    def merge(self):
        """This function processes multiple percolator runs from fractions and
        filters out the best scoring peptides. It writes a single fraction with
        those peptides and ALL psms from all fractions. 
        
        Namespace and static xml come from first percolator file. 
        Make sure fractions are from same percolator run."""
        self.prepare_merge()
        uniquepeps = filtering.filter_unique_peptides(self.allpeps, self.score,
                                                        self.ns)
        self.features = {'psm': self.allpsms, 'peptide': uniquepeps}

