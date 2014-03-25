import os

from app.lookups import quant as lookups
from app.readers import openmz as readers
from app.preparation import openmz as preparation
from app.writers import tsv as writers


class BaseDriver(object):
    def __init__(self, **kwargs):
        self.outdir = kwargs['outdir']

    def create_outfilepath(self, fn, suffix=None):
        basefn, ext = os.path.splitext(os.path.basename(fn))
        outfn = basefn + suffix + ext
        return os.path.join(self.outdir, outfn)


class TSVQuantDriver(BaseDriver):
    def __init__(self, **kwargs):
        super(TSVQuantDriver, self).__init__(**kwargs)
        self.spectrafn = kwargs.get('spectra', None)
        self.quantfn = kwargs.get('quants', None)
        self.tsvfn = kwargs.get('tsv', None)

    def run(self):
        quantdb = self.create_quantlookup()
        self.generate_psms(quantdb)
        self.write()

    def create_quantlookup(self):
        """Creates sqlite file containing quantification data and
        spectra file/scan number data for looking up quant data.
        Returns sqlite file name.
        """
        spectra = readers.mzml_generator(self.spectrafn)
        consensus_quants = readers.quant_generator(self.quantfn)
        return lookups.create_quant_lookup(spectra, consensus_quants)

    def generate_psms(self, quantdb):
        self.psms = preparation.generate_psms_quanted(quantdb, self.tsvfn)

    def write(self):
        outfn = self.create_outfilepath(self.tsvm, 'quants')
        writers.write_quantpsm_tsv(self.header, self.psms, outfn)
