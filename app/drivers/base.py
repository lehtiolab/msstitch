import os
import sys

from app.lookups import base as lookups
from app.drivers.options import shared_options


class BaseDriver(object):
    def __init__(self):
        self.lookupfn = None
        self.shared_options = shared_options
        self.options = ['fn', 'outdir']

    def set_lookup(self):
        if self.lookupfn is not None and hasattr(self, 'lookuptype'):
            self.lookup = lookups.get_lookup(self.lookupfn, self.lookuptype)
        else:
            self.lookup = None

    def set_options(self, options=None):
        if options is not None:
            self.options.extend(options)
        self.options = self.get_parser_options(self.options)

    def get_commandhelp(self):
        return self.commandhelp

    def get_options(self):
        return self.options

    def get_parser_options(self, names):
        """Given a list of option names, this returns a list of dicts
        defined in all_options and self.shared_options. These can then
        be used to populate the argparser with"""
        options = []
        for name in names:
            try:
                option = self.parser_options[name]
            except KeyError:
                option = self.shared_options[name]
            options.append(option)
        return options

    def parse_input(self, **kwargs):
        for option in self.options:
            opt_argkey = option['driverattr']
            if 'type' in option and option['type'] == 'pick':
                try:
                    assert kwargs.get(opt_argkey) in option['picks']
                except AssertionError:
                    print('Option {} should be one of [{}]'.format(
                        option['clarg'], ','.join(option['picks'])))
                    sys.exit(1)
            setattr(self, opt_argkey, kwargs.get(opt_argkey))

    def start(self, **kwargs):
        self.parse_input(**kwargs)
        self.set_lookup()
        self.run()

    def finish(self):
        """Cleans up after work"""
        pass

    def create_outfilepath(self, fn, suffix=None):
        basefn = os.path.basename(fn)
        outfn = basefn + suffix
        return os.path.join(self.outdir, outfn)

    def create_multi_outfile_basepath(self, fn, suffix=None):
        """Returns a basepath to dynamically create outfiles by writer
        modules. Basepath includes a formatting string"""
        return self.create_outfilepath(fn + '_{0}', suffix)

    def get_column_header_for_number(self, column_var_names, header=False):
        """This function subtracts 1 from inputted column number to comply
        with programmers counting (i.e. from 0, not from 1). For TSV data."""
        if not header:
            header = self.oldheader
        for col in column_var_names:
            value = getattr(self, col)
            if value is None:
                continue
            setattr(self, col, header[int(value) - 1])

    def get_commandname(self):
        return self.command
