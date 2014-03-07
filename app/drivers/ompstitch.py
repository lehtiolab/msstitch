import os

from app.lookups import quant as lookups
from app.readers import ompstitch as readers
from app.preparation import ompstitch as preparation
from app.writers import ompstitch as writers


class BaseDriver(object):
    def __init__(self, **kwargs):
        self.outdir = kwargs['outdir']

    def create_outfilepath(self, fn, suffix=None):
        basefn = os.path.basename(fn)
        outfn = basefn + suffix
        return os.path.join(self.outdir, outfn)


class QuantTSVDriver(BaseDriver):
    def __init__(self, **kwargs):
        super(QuantTSVDriver, self).__init__(**kwargs)
        self.spectrafn = kwargs.get('spectra', None)
        self.quantfn = kwargs.get('quants', None)
        self.idfn = kwargs.get('ids', None)
        if None in [self.spectrafn, self.quantfn, self.idfn]:
            raise Exception(
                'Need to specify spectra file, quant file and PSM file')

    def run(self):
        self.get_psms()
        self.write()

    def get_psms(self):
        """Creates a lookup db with scan nrs and quants. Then looks up psms
        in the db via their scan nr and outputs generator of psms with
        quant info."""
        spectra = readers.mzml_generator(self.spectrafn)
        consensus_quants = readers.quant_generator(self.quantfn)
        quantdbfn = lookups.create_quant_lookup(spectra, consensus_quants)
        psms = readers.psm_generator(self.idfn)
        self.psms = preparation.generate_psms_quanted(quantdbfn, psms)

    def write(self):
        writers.write_quantpsm_tsv(self.psms)
