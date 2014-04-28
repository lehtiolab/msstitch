import os

from app.readers import pycolator as readers
from app.writers import pycolator as writers


class BaseDriver(object):
    def __init__(self, **kwargs):
        self.fn = kwargs['infile']
        self.outdir = kwargs['outdir']

    def create_outfilepath(self, fn, suffix=None):
        basefn = os.path.basename(fn)
        outfn = basefn + suffix
        return os.path.join(self.outdir, outfn)


class PycolatorDriver(BaseDriver):
    """Driver for pycolator functions"""
    def prepare_percolator_output(self, fn):
        """Returns namespace and static xml from percolator output file"""
        ns = readers.get_namespace(fn)
        static = readers.get_percolator_static_xml(fn, ns)
        return ns, static

    def get_all_peptides(self):
        return readers.generate_peptides(self.fn, self.ns)

    def get_all_psms(self):
        return readers.generate_psms(self.fn, self.ns)

    def prepare(self, fn):
        self.ns, self.static_xml = self.prepare_percolator_output(fn)
        self.allpeps = self.get_all_peptides()
        self.allpsms = self.get_all_psms()

    def run(self):
        self.prepare(self.fn)
        self.set_features()
        self.write()

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_percolator_xml(self.static_xml, self.features, outfn)
