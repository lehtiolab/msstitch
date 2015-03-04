import shutil
from app.readers import spectra as spectrareader
from app.readers import openms as openmsreader
from app.lookups import quant as lookups
from app.drivers.base import BaseDriver


class LookupDriver(BaseDriver):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.lookupfn = kwargs.get('lookup', None)
        if self.lookupfn is not None:
            # FIXME make this general
            self.lookup = lookups.get_lookup(self.lookupfn, self.lookuptype)
        else:
            self.lookupfn = 'msstitcher_lookup.sqlite'
            self.lookup = lookups.initiate_quant_lookup()

    def run(self):
        self.create_lookup()

    def write_move(self):
        """Moves outfile from workdir to destination, used from different
        lookup creating commands"""
        outfn = self.create_outfilepath(self.lookupfn, self.outsuffix)
        shutil.move(self.lookup.get_fn(), outfn)


class QuantLookupDriver(LookupDriver):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.spectrafns = kwargs.get('spectra', None)  # not for all lookups


class SpectraLookupDriver(QuantLookupDriver):
    outsuffix = '_spectralookup.sqlite'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.spectrafns = self.fn

    def create_lookup(self):
        fn_spectra = spectrareader.mzmlfn_spectra_generator(self.spectrafns)
        lookups.create_spectra_lookup(self.lookup, fn_spectra)


class IsobaricQuantLookupDriver(QuantLookupDriver):
    outsuffix = '_isobquantlookup.sqlite'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.consensusfns = self.fn

    def create_lookup(self):
        # FIXME here get quantmap for channel labels as dict {'0': '113', etc}
        # then pass this dict to create_isobaric_quant_lookup
        quantmap = openmsreader.get_quantmap(self.consensusfns)
        mzmlfn_consxml = openmsreader.mzmlfn_cons_el_generator(self.spectrafns,
                                                               self.consensusfns)
        lookups.create_isobaric_quant_lookup(self.lookup, mzmlfn_consxml,
                                             quantmap),


class PrecursorQuantLookupDriver(QuantLookupDriver):
    outsuffix = '_ms1quantlookup.sqlite'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.precursorfns = self.fn

    def create_lookup(self):
        specfn_feats = openmsreader.mzmlfn_feature_generator(self.spectrafns,
                                                             self.precursorfns)
        lookups.create_precursor_quant_lookup(self.lookup, specfn_feats)
