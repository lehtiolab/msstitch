import shutil
from app.readers import spectra as spectrareader
from app.readers import openms as openmsreader
from app.lookups import quant as lookups
from app.drivers import BaseDriver


class LookupDriver(BaseDriver):
    def __init__(self, **kwargs):
        lookupfn = kwargs.get('lookup', None)
        if lookupfn is not None:
            # FIXME make this general
            self.lookup = lookups.get_quant_lookup(lookupfn)
        else:
            self.lookup = lookups.initiate_quant_lookup(self.workdir,
                                                        foreign_keys=True)

    def run(self):
        self.create_lookup()
        self.write_move()
        self.finish()

    def write_move(self):
        """Moves outfile from workdir to destination, used from different
        lookup creating commands"""
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        shutil.move(self.lookup.get_fn(), outfn)


class QuantLookupDriver(LookupDriver):
    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.spectrafns = kwargs.get('spectra', None)  # not for all lookups


class SpectraLookupDriver(QuantLookupDriver):
    outsuffix = 'spectralookup.sqlite'

    def create_lookup(self):
        fn_spectra = spectrareader.mzmlfn_spectra_generator(self.spectrafns)
        lookups.create_spectra_lookup(self.lookup, fn_spectra)


class IsobaricQuantLookupDriver(QuantLookupDriver):
    outsuffix = 'isobquantlookup.sqlite'

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.consensusfns = kwargs.get('isob_quants', None)

    def create_lookup(self):
        # FIXME here get quantmap for channel labels as dict {'0': '113', etc}
        # then pass this dict to create_isobaric_quant_lookup
        mzmlfn_consxml = openmsreader.mzmlfn_cons_el_generator(self.spectrafns,
                                                               self.quantfns)
        lookups.create_isobaric_quant_lookup(self.quantdb, mzmlfn_consxml),


class PrecursorQuantLookupDriver(QuantLookupDriver):
    outsuffix = '_ms1quantlookup.sqlite'

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.precursorfns = kwargs.get('precur_quants', None)

    def create_lookup(self):
        specfn_feats = openmsreader.mzmlfn_feature_generator(self.spectrafns,
                                                             self.precursorfns)
        lookups.create_precursor_quant_lookup(self.quantdb, specfn_feats)
