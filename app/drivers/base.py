import os
import shutil
from tempfile import mkdtemp

from app.lookups import base as lookups


class BaseDriver(object):
    def __init__(self, **kwargs):
        self.fn = kwargs['infile']
        self.outdir = kwargs['outdir']
        lookupfn = kwargs.get('lookup', None)
        if lookupfn is not None and hasattr(self, 'lookuptype'):
            self.lookup = lookups.get_lookup(lookupfn, self.lookuptype)
        else:
            self.lookup = None
        # self.workdir = self.set_workdir(kwargs.get('workdir', os.getcwd()))

    def finish(self):
        """Cleans up after work"""
        pass
        #self.clean_workdir()

    def set_workdir(self, parentworkdir):
        """Creates temporary workdir and returns absolute path to it"""
        # FIXME DEPRECATE no more workdirs
        return mkdtemp(suffix='_msstitcher', dir=parentworkdir)

    def clean_workdir(self):
        """Removes workdir from filesystem"""
        # FIXME DEPRECATE no more workdirs
        shutil.rmtree(self.workdir)

    def create_outfilepath(self, fn, suffix=None):
        basefn = os.path.basename(fn)
        outfn = basefn + suffix
        return os.path.join(self.outdir, outfn)

    def create_multi_outfile_basepath(self, fn, suffix=None):
        """Returns a basepath to dynamically create outfiles by writer
        modules. Basepath includes a formatting string"""
        return self.create_outfilepath(fn + '_{0}', suffix)
