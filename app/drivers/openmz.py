from app.drivers.openmz import base
from app.lookups import quant as lookups
from app.readers import openmz as readers
from app.preparation import openmz as preparation
from app.writers import tsv as writers


class TSVQuantDriver(base.OpenMzDriver):
    def __init__(self, **kwargs):
        super(TSVQuantDriver, self).__init__(**kwargs)
        self.spectrafns = kwargs.get('spectra', None)
        self.quantfn = kwargs.get('infile', None)
        self.tsvfn = kwargs.get('mzidtsv', None)

    def run(self):
        quantdb = self.create_quantlookup()
        self.generate_psms(quantdb)
        self.write()

    def create_quantlookup(self):
        """Creates sqlite file containing quantification data and
        spectra file/scan number data for looking up quant data.
        Returns sqlite file name.
        """
        fn_spectra = readers.mzml_generator(self.spectrafns)
        consensus_quants = readers.quant_generator(self.quantfn)
        return lookups.create_quant_lookup(fn_spectra, consensus_quants)

    def generate_tsv_content(self, quantdb):
        """Creates iterator to write to new tsv. Contains input tsv
        lines plus quant data for these."""
        quantheader = preparation.get_quant_header(quantdb)
        self.header = preparation.create_tsv_header_quant(self.tsvfn,
                                                          quantheader)
        self.psms = preparation.generate_psms_quanted(quantdb,
                                                      self.tsvfn,
                                                      quantheader)

    def write(self):
        outfn = self.create_outfilepath(self.tsvfn, self.outsuffix)
        writers.write_quantpsm_tsv(self.header, self.psms, outfn)
