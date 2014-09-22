from app.preparation import mzidplus as prep
from app.drivers.mzidtsv import MzidTSVDriver


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
