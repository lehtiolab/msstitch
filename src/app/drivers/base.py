import os
import sys

from app.lookups import base as lookups
from app.drivers.options import (shared_options, lookup_options, psmtable_options)

from app.readers import tsv as tsvreader
from app.readers import percolator as percoreaders
from app.readers import xml

from app.writers import tsv as tsvwriter
from app.writers import percolator as percowriters


class BaseDriver(object):
    tabletypes = []

    def __init__(self):
        self.lookupfn = None
        self.infiletype = ''

    def set_lookup(self):
        if self.lookupfn is not None and hasattr(self, 'lookuptype'):
            self.lookup = lookups.get_lookup(self.lookupfn, self.lookuptype)
        else:
            self.lookup = None

    def initialize_lookup(self, outfile=None):
        if self.lookup is None:
            # FIXME MUST be a set or mzml lookup? here is place to assert
            # correct lookuptype!
            if outfile is None and self.outfile is None:
                self.outfile = os.path.join(self.outdir,
                                            'mslookup_db.sqlite')
                lookupfile = self.outfile
            elif outfile is not None:
                lookupfile = outfile
            elif self.outfile is not None:
                lookupfile = self.outfile
            self.lookup = lookups.create_new_lookup(lookupfile,
                                                    self.lookuptype)
        self.lookup.add_tables(self.tabletypes)

    def set_options(self):
        self.options = self.define_options(['fn', 'outdir', 'outfile'], {})
        #self.options['-i']['help'] = self.options['-i']['help'].format(
        #@    self.infiletype)

    def get_commandhelp(self):
        return self.commandhelp

    def define_options(self, names, parser_options=None):
        """Given a list of option names, this returns a list of dicts
        defined in all_options and self.shared_options. These can then
        be used to populate the argparser with"""
        if parser_options is None:
            parser_options = {}
        options = {}
        for name in names:
            try:
                option = {k: v for k, v in parser_options[name].items()}
            except KeyError:
                option = {k: v for k, v in shared_options[name].items()}
            try:
                options.update({name: option})
                #options.update({option['clarg']: option})
            except TypeError:
                #options.update({option['clarg'][0]: option})
                options.update({name: option})
        return options

    def parse_input(self, **kwargs):
        # Set option values
        for option in self.options.values():
            opt_argkey = option['driverattr']
            opt_val = kwargs.get(opt_argkey)
            setattr(self, opt_argkey, opt_val)
        # Parse conditional required options
        required = {}
        for optname, option in self.options.items():
            key = option['driverattr']
            conds = option.get('conditional_required', [])
            for cond in conds:
                if getattr(self, self.options[cond]['driverattr']) and not getattr(self, key):
                    try:
                        required[cond].append(optname) 
                    except KeyError:
                        required[cond] = [optname]
        for cond, keys in required.items():
            print('When using {}, you must also specify {}'.format(
                self.options[cond]['clarg'], ', '.join([self.options[k]['clarg'] for k in keys])))
            sys.exit(1)
        # Prepare for output
        if self.outdir is None:
            self.outdir = os.getcwd()
        if self.outfile is not None:
            if not os.path.isabs(self.outfile):
                self.outfile = os.path.join(self.outdir, self.outfile)

    def start(self, **kwargs):
        self.parse_input(**kwargs)
        self.set_lookup()
        self.run()

    def create_outfilepath(self, fn, suffix=None):
        if self.outfile is None:
            basefn = os.path.basename(fn)
            outfn = os.path.join(self.outdir, basefn + suffix)
        else:
            outfn = self.outfile
        return outfn

    def number_to_headerfield(self, columnr, header):
        return header[int(columnr) - 1]

    def get_column_header_for_number(self, column_var_names, header=False):
        """This function subtracts 1 from inputted column number to comply
        with programmers counting (i.e. from 0, not from 1). For TSV data."""
        if not header:
            header = self.oldheader
        for col in column_var_names:
            value = getattr(self, col)
            if not value or value is None:
                continue
            setattr(self, col, self.number_to_headerfield(value, header))

    def get_commandname(self):
        return self.command

    def get_fastadelim_genefield(self, delimtype, genefield):
        if delimtype is not None:
            fastadelim = {'tab': '\t',
                          'semicolon': ';', 'pipe': '|'}[delimtype]
        else:
            fastadelim = None
        if genefield is not None:
            genefield -= 1
        return fastadelim, genefield

    def prepare(self):
        pass

    def run(self):
        self.prepare()
        self.set_features()
        self.write()


class LookupDriver(BaseDriver):
    outfile, outdir = None, None
    tabletypes = []

    def __init__(self):
        super().__init__()
        self.parser_options = lookup_options

    def set_options(self):
        super().set_options()
        del(self.options['fn'])
        del(self.options['outfile'])
        del(self.options['outdir'])
        self.options.update(self.define_options(['lookupfn'],
                                                lookup_options))

    def run(self):
        self.initialize_lookup()
        self.create_lookup()


class PSMDriver(BaseDriver):
    def __init__(self):
        super().__init__()
        self.infiletype = 'TSV PSM table (MSGF+)'
        self.parser_options = psmtable_options

    def prepare(self):
        if type(self.fn) == list:
            self.first_infile = self.fn[0]
        else:
            self.first_infile = self.fn
        self.oldheader = tsvreader.get_tsv_header(self.first_infile)
        self.oldpsms = tsvreader.generate_split_tsv_lines(self.fn, self.oldheader)

    def write(self):
        outfn = self.create_outfilepath(self.first_infile, self.outsuffix)
        tsvwriter.write_tsv(self.header, self.psms, outfn)


class PercolatorDriver(BaseDriver):
    """Driver for percolator functions"""
    def __init__(self):
        super().__init__()
        self.infiletype = 'percolator out XML'

    def prepare_percolator_output(self, fn):
        """Returns namespace and static xml from percolator output file"""
        ns = xml.get_namespace(fn)
        static = percoreaders.get_percolator_static_xml(fn, ns)
        return ns, static

    def get_all_peptides(self):
        return percoreaders.generate_peptides(self.fn, self.ns)

    def get_all_psms(self):
        return percoreaders.generate_psms(self.fn, self.ns)

    def get_all_psms_strings(self):
        return percoreaders.generate_psms_multiple_fractions_strings([self.fn],
                                                                self.ns)

    def get_all_peptides_strings(self):
        return percoreaders.generate_peptides_multiple_fractions_strings([self.fn],
                                                                    self.ns)

    def prepare(self):
        self.ns, self.static_xml = self.prepare_percolator_output(self.fn)
        self.allpeps = self.get_all_peptides()
        self.allpsms = self.get_all_psms()

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        percowriters.write_percolator_xml(self.static_xml, self.features, outfn)


class PepProttableDriver(BaseDriver):

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        tsvwriter.write_table_with_na(self.header, self.features, outfn)
