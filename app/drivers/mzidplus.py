from app.drivers import base
from app.preparation import mzidplus as prep
from app.preparation import quant as quantprep
from app.writers import mzidplus as writers
from app.readers import tsv as tsvreader
from app.readers import spectra as spectrareader
from app.readers import openms as openmsreader
from app.lookups import quant as lookups


class MzidPlusDriver(base.BaseDriver):
    def run(self):
        self.get_psms()
        self.write()

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_mzid_tsv(self.header, self.psms, outfn)


class MzidTSVConcatenateDriver(MzidPlusDriver):
    """Concatenates TSVs"""
    outsuffix = '_concat.tsv'

    def __init__(self, **kwargs):
        super(MzidTSVConcatenateDriver, self).__init__(**kwargs)
        self.allinfiles = [self.fn]
        self.allinfiles.extend(kwargs.get('multifile_input', None))

    def get_psms(self):
        self.header = tsvreader.get_tsv_header(self.fn)
        self.psms = prep.merge_mzidtsvs(self.allinfiles, self.header)


class MzidPercoTSVDriver(MzidPlusDriver):
    """
    Adds percolator data from mzid file to table.
    """
    outsuffix = '_percolated.tsv'

    def __init__(self, **kwargs):
        super(MzidPercoTSVDriver, self).__init__(**kwargs)
        self.idfn = kwargs.get('mzid', None)
        self.multipsm_per_scan = kwargs.get('allpsms', False)
        assert self.idfn is not None

    def get_psms(self):
        if self.multipsm_per_scan is True:
            # FIXME not supported yet
            # Create mzid PSM/sequence sqlite (fn, scan, rank, sequence)
            pass
        else:
            seqlookup = None

        oldheader = tsvreader.get_tsv_header(self.fn)
        self.header = prep.get_header_with_percolator(oldheader,
                                                      self.multipsm_per_scan)
        self.psms = prep.add_percolator_to_mzidtsv(self.idfn,
                                                   self.fn,
                                                   self.multipsm_per_scan,
                                                   oldheader,
                                                   self.header,
                                                   seqlookup)


class TSVQuantDriver(MzidPlusDriver):
    outsuffix = '_quant.tsv'

    def __init__(self, **kwargs):
        super(TSVQuantDriver, self).__init__(**kwargs)
        self.spectrafns = kwargs.get('spectra', None)
        self.quantfns = kwargs.get('quants', None)

    def get_psms(self):
        """Creates iterator to write to new tsv. Contains input tsv
        lines plus quant data for these."""
        quantdb = self.create_quantlookup()
        oldheader = tsvreader.get_tsv_header(self.fn)
        self.header, qheader = quantprep.get_full_and_quant_headers(oldheader,
                                                                    quantdb)
        self.psms = quantprep.generate_psms_quanted(quantdb, self.fn,
                                                    qheader, oldheader)

    def create_quantlookup(self):
        """Creates sqlite file containing quantification data and
        spectra file/scan number data for looking up quant data.
        Returns sqlite file name.
        """
        fn_spectra = spectrareader.mzml_generator(self.spectrafns)
        consensus_quants = openmsreader.quant_generator(self.quantfns)
        return lookups.create_quant_lookup(fn_spectra, consensus_quants)
