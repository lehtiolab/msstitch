from app.drivers import base
from app.readers import pycolator as readers
from app.readers import xml
from app.writers import pycolator as writers


class PycolatorDriver(base.BaseDriver):
    """Driver for pycolator functions"""
    def __init__(self):
        super().__init__()
        self.infiletype = 'percolator out XML'

    def prepare_percolator_output(self, fn):
        """Returns namespace and static xml from percolator output file"""
        ns = xml.get_namespace(fn)
        static = readers.get_percolator_static_xml(fn, ns)
        return ns, static

    def get_all_peptides(self):
        return readers.generate_peptides(self.fn, self.ns)

    def get_all_psms(self):
        return readers.generate_psms(self.fn, self.ns)

    def get_all_psms_strings(self):
        return readers.generate_psms_multiple_fractions_strings([self.fn],
                                                                self.ns)

    def get_all_peptides_strings(self):
        return readers.generate_peptides_multiple_fractions_strings([self.fn],
                                                                    self.ns)

    def prepare(self):
        self.ns, self.static_xml = self.prepare_percolator_output(self.fn)
        self.allpeps = self.get_all_peptides()
        self.allpsms = self.get_all_psms()

    def run(self):
        self.prepare()
        self.set_features()
        self.write()
        self.finish()

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_percolator_xml(self.static_xml, self.features, outfn)
