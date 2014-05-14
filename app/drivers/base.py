import os


class BaseDriver(object):
    def __init__(self, **kwargs):
        self.fn = kwargs['infile']
        self.outdir = kwargs['outdir']

    def create_outfilepath(self, fn, suffix=None):
        basefn = os.path.basename(fn)
        outfn = basefn + suffix
        return os.path.join(self.outdir, outfn)
