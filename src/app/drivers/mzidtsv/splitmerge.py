import os

from app.actions.mzidtsv import splitmerge as prep
from app.drivers.mzidtsv import MzidTSVDriver
from app.writers import mzidtsv as writers
from app.drivers.options import mzidtsv_options


class MzidTSVConcatenateDriver(MzidTSVDriver):
    """Concatenates TSVs"""
    outsuffix = '_concat.tsv'
    command = 'merge'
    commandhelp = ('Merges multiple TSV tables of MSGF+ output.'
                   'Make sure headers are same in all files.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['multifiles'],
                                                mzidtsv_options))

    def get_psms(self):
        self.header = self.oldheader
        self.psms = prep.merge_mzidtsvs(self.fn, self.oldheader)


class MzidTSVSplitDriver(MzidTSVDriver):
    """Splits MSGF PSM table on contents of certain column. Each
    row in file is piped to an output file. Which output files
    the row is written to depends on the contents of the selected
    column"""
    outsuffix = '_split.tsv'
    command = 'split'
    commandhelp = 'Splits an MSGF TSV PSM table into multiple new tables'

    def set_options(self):
        super().set_options()
        options = self.define_options(['bioset', 'splitcol'],
                                      mzidtsv_options)
        self.options.update(options)

    def get_psms(self):
        self.header = self.oldheader[:]
        self.psms = prep.generate_psms_split(self.fn, self.oldheader,
                                             self.bioset, self.splitcol)

    def write(self):
        base_outfile = os.path.join(self.outdir, '{}.tsv')
        writers.write_multi_mzidtsv(self.header, self.oldheader, self.psms,
                                    base_outfile)
