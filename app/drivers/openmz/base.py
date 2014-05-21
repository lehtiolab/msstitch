from app.drivers import base
from app.readers import basereader
from app.writers import tsv as writers


class OpenMzDriver(base.BaseDriver):
    """Driver for OpenMz functions"""

    def run(self):
        self.prepare()
        self.set_features()
        self.write()

    def prepare(self):
        self.ns = basereader.get_namespace(self.fn)

    def write(self):
        outfn = self.create_outfilepath(self.tsvfn, self.outsuffix)
        writers.write_quantpsm_tsv(self.header, self.psms, outfn)
