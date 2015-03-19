from app.actions.mzidtsv import splitmerge as prep
from app.drivers.mzidtsv import MzidTSVDriver
from app.writers import mzidtsv as writers


class MzidTSVConcatenateDriver(MzidTSVDriver):
    """Concatenates TSVs"""
    outsuffix = '_concat.tsv'

    def __init__(self, **kwargs):
        super(MzidTSVConcatenateDriver, self).__init__(**kwargs)
        self.allinfiles = [self.fn]
        self.allinfiles.extend(kwargs.get('multifile_input', None))

    def get_psms(self):
        self.header = self.oldheader
        self.psms = prep.merge_mzidtsvs(self.allinfiles, self.oldheader)


class MzidTSVSplitDriver(MzidTSVDriver):
    """Splits MSGF PSM table on contents of certain column. Each
    row in file is piped to an output file. Which output files
    the row is written to depends on the contents of the selected
    column"""
    outsuffix = '_split.tsv'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.bioset = kwargs.get('bioset', None)
        self.splitcol = kwargs.get('splitcol', None)
        self.renamecols = kwargs.get('renamecols', None)
        self.renamecolpattern = kwargs.get('renamecolpattern', None)

    def get_psms(self):
        self.header = prep.get_splitheader(self.oldheader, self.renamecols,
                                           self.renamecolpattern)
        self.psms = prep.generate_psms_split(self.fn, self.oldheader,
                                             self.bioset, self.splitcol)

    def write(self):
        base_outfile = self.create_multi_outfile_basepath(self.fn, self.outsuffix)
        writers.write_multi_mzidtsv(self.header, self.oldheader, self.psms, base_outfile)
