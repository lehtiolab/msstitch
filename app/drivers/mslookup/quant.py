from app.drivers.mslookup import base
from app.readers import openms as openmsreader
from app.readers import tsv as tsvreader
from app.actions.mslookup import quant as lookups
from app.drivers.options import mslookup_options


class QuantLookupDriver(base.LookupDriver):
    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['spectrafns', 'multifiles'],
                                                mslookup_options))


class IsobaricQuantLookupDriver(QuantLookupDriver):
    lookuptype = 'isoquant'
    command = 'isoquant'
    commandhelp = ('Create lookup of isobaric quant data in OpenMS '
                   'consensusXML format. Use requires --spectra, --dbfile '
                   'with an sqlite lookup of spectra, and passing multiple '
                   'consensusXML files to -i.')

    def __init__(self):
        super().__init__()
        self.infiletype = 'consensusXML'

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.cons_fns = self.fn

    def create_lookup(self):
        quantmap = openmsreader.get_quantmap(self.cons_fns[0])
        mzmlfn_consxml = openmsreader.mzmlfn_cons_el_generator(self.spectrafns,
                                                               self.cons_fns)
        lookups.create_isobaric_quant_lookup(self.lookup, mzmlfn_consxml,
                                             quantmap),


class PrecursorQuantLookupDriver(QuantLookupDriver):
    lookuptype = 'ms1quant'
    command = 'ms1quant'
    commandhelp = ('Creates lookup of precursor quant data in OpenMS '
                   'featureXML, or Kronik output formats. '
                   'Use requires --spectra, --dbfile with an sqlite lookup of '
                   'spectra, --quanttype to determine quant output, --mztol, '
                   '--rttol, --mztoltype for tolerance specification, and '
                   'passing a featureXML or kronik file to -i.')

    def __init__(self):
        super().__init__()
        self.infiletype = 'featureXML or Kronik'

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['quantfiletype', 'rttol',
                                                 'mztol', 'mztoltype'],
                                                mslookup_options))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.precursorfns = self.fn

    def create_lookup(self):
        if self.quantfiletype == 'openms':
            specfn_feats = openmsreader.mzmlfn_feature_generator(
                self.spectrafns, self.precursorfns)
        elif self.quantfiletype == 'kronik':
            specfn_feats = tsvreader.mzmlfn_kronikfeature_generator(
                self.spectrafns, self.precursorfns)
        lookups.create_precursor_quant_lookup(self.lookup, specfn_feats,
                                              self.quantfiletype,
                                              self.rt_tol,
                                              self.mz_tol,
                                              self.mz_toltype)
