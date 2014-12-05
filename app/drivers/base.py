import os
import shutil
from tempfile import mkdtemp


class BaseDriver(object):
    def __init__(self, **kwargs):
        self.fn = kwargs['infile']
        self.outdir = kwargs['outdir']
        self.workdir = self.set_workdir(kwargs.get('workdir', os.getcwd()))

    def finish(self):
        """Cleans up after work"""
        self.clean_workdir()

    def set_workdir(self, parentworkdir):
        """Creates temporary workdir and returns absolute path to it"""
        return mkdtemp(suffix='_msstitcher', dir=parentworkdir)

    def clean_workdir(self):
        """Removes workdir from filesystem"""
        shutil.rmtree(self.workdir)

    def create_outfilepath(self, fn, suffix=None):
        basefn = os.path.basename(fn)
        outfn = basefn + suffix
        return os.path.join(self.outdir, outfn)

    def copy_db_to_workdir(self):
        lookupfn = os.path.basename(self.lookup)
        copied_lookup = os.path.join(self.workdir, lookupfn)
        shutil.copy(self.lookup, copied_lookup)
        self.lookup = copied_lookup
