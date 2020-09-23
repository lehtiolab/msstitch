from app.drivers import base
from app.drivers.options import lookup_options, sequence_options

from app.readers import spectra as spectrareader
from app.readers import openms as openmsreader
from app.readers import tsv as tsvreader

from app.actions.lookups import spectra as spectralookup
from app.actions.lookups import quant as quantlookups 
from app.actions.lookups import sequence as seqlookups


# FIXME infn/outfn is special and should be speicifed on driver, or on baselookup?
# Or is it only used in storeseq?


class SpectraLookupDriver(base.LookupDriver):
    lookuptype = 'spectra'
    command = 'storespectra'
    commandhelp = ('Create lookup of spectra in mzML format. '
            'Requires passing mzML files to --spectra. '
            'Biological set names for each file should be specified using '
            '--setnames')

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        if self.setnames is None:
            assert self.lookup is not None, ('Must specify lookup '
                                             'if setnames have not '
                                             'been provided')
        else:
            self.setnames = [x.replace('"', '') for x in self.setnames]

    def set_options(self):
        super().set_options()
        self.options['lookupfn'].update({'required': False, 'default': None})
        self.options.update(self.define_options(['outfile', 'spectrafns', 'setnames'],
                                                lookup_options))

    def create_lookup(self):
        spectralookup.create_bioset_lookup(self.lookup, self.spectrafns,
                                          self.setnames)
        fn_spectra = spectrareader.mzmlfn_ms2_spectra_generator(
            self.spectrafns)
        spectralookup.create_spectra_lookup(self.lookup, fn_spectra)


class QuantLookupDriver(base.LookupDriver):
    command = 'storequant'
    lookuptype = 'specquant'

    commandhelp = ('Create lookup of isobaric quant data in OpenMS ' # FIXME
                   'consensusXML format. Use requires matched files passed to '
                   '--spectra, --ms1, --isobaric parameters, and --dbfile '
                   'with an sqlite lookup of already stored spectra.'
                   'Creates lookup of precursor quant data in OpenMS '
                   'featureXML, or Kronik output formats. '
                   'Use requires --spectra, --dbfile with an sqlite lookup of '
                   'spectra, --quanttype to determine quant output, --mztol, '
                   '--rttol, --mztoltype for tolerance specification, and '
                   'passing a featureXML or kronik file to -i.')



    def set_options(self):
        super().set_options()
        self.options.update(self.define_options([
            'spectrafns', 'kronik', 'dinosaur', 'isobaric',
            'sum_or_apex', 'rttol', 'mztol', 'mztoltype'], lookup_options))

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        if getattr(self, 'isobaricfns') and len(self.isobaricfns) > 0:
            self.tabletypes.append('isobaric')
        else:
            self.isobaricfns = False
        self.ms1type = False
        if getattr(self, 'dinosaurfns') and len(self.dinosaurfns) > 0:
            self.tabletypes.append('ms1')
            self.ms1type = 'dinosaur'
            self.ms1fns = self.dinosaurfns[:]
        elif getattr(self, 'kronikfns') and len(self.kronikfns) > 0:
            self.tabletypes.append('ms1')
            self.ms1type = 'kronik'
            self.ms1fns = self.kronikfns[:]

    def create_lookup(self):
        if self.isobaricfns:
            quantmap = openmsreader.get_quantmap(self.isobaricfns[0])
            mzmlfn_consxml = openmsreader.mzmlfn_cons_el_generator(self.spectrafns,
                    self.isobaricfns)
            quantlookups.create_isobaric_quant_lookup(self.lookup, mzmlfn_consxml, 
                    quantmap),
        if self.ms1type:
            ms1feats = tsvreader.mzmlfn_tsvfeature_generator(self.spectrafns, self.ms1fns)
            quantlookups.create_precursor_quant_lookup(self.lookup, ms1feats,
                    self.sum_or_apex, self.ms1type, self.rt_tol, self.mz_tol, self.mz_toltype)


class SequenceLookupDriver(base.LookupDriver):
    """Creates an SQLite lookup DB from a FASTA file. Sequences are
    trypsinized and stored. It's possible to store sequences reversed
    for N-terminal falloff indexing, and it can be specified to cut
    tryptic before proline.
    """
    lookuptype = 'searchspace'
    command = 'storeseq'
    commandhelp = """Create a lookup DB from a FASTA file. Sequences are
    either stored as 1-aa-shifted overlapping peptides for full-length-protein 
    matching, or as trypsinized peptides."""

    def __init__(self):
        super().__init__()
        self.infiletype = 'FASTA'

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['fn', 'outfile']))
        self.options['lookupfn'].update({'required': False, 'default': None})
        self.options.update(self.define_options(['fullprotein', 'falloff',
            'proline', 'minlength'], lookup_options))
        self.options.update(self.define_options(['trypsinize', 'miss_cleavage'],
            sequence_options))

    def create_lookup(self):
        if self.fullprotein:
            if not self.minlength:
                self.minlength = 7
                # Default here since it is a param shared with many others 
                # where it can be 0. TODO create own param for this driver
                print('--minlen not specified, default to 7')
            print('Creating full-length protein lookup, minimum matching '
                    'peptide length is {}'.format(self.minlength))
            if self.proline or self.falloff or not self.trypsinize or self.miss_cleavage:
                print('Ignoring other options for tryptic lookup building')
            seqlookups.create_searchspace_wholeproteins(self.lookup, self.fn,
                                                         self.minlength)
        else:
            seqlookups.create_searchspace(self.lookup, self.fn, self.minlength, self.proline,
                                           self.falloff, self.trypsinize, self.miss_cleavage)
