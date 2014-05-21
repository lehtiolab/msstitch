from app.drivers.openmz import base
from app.lookups import quant as lookups
from app.readers import openmz as readers
from app.preparation import openmz as preparation


class TSVQuantDriver(base.OpenMzDriver):
    outsuffix = '_quant.tsv'

    def __init__(self, **kwargs):
        super(TSVQuantDriver, self).__init__(**kwargs)
        self.spectrafns = kwargs.get('spectra', None)
        self.quantfn = kwargs.get('infile', None)
        self.tsvfn = kwargs.get('mzidtsv', None)

    def set_features(self):
        quantdb = self.create_quantlookup()
        self.generate_psms(quantdb)

    def create_quantlookup(self):
        """Creates sqlite file containing quantification data and
        spectra file/scan number data for looking up quant data.
        Returns sqlite file name.
        """
        fn_spectra = readers.mzml_generator(self.spectrafns)
        consensus_quants = readers.quant_generator(self.quantfn, self.ns)
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
