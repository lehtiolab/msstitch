from app.preparation.mzidtsv import quant as prep
from app.drivers.mzidtsv import MzidTSVDriver


class TSVQuantDriver(MzidTSVDriver):
    outsuffix = '_quant.tsv'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.quantdb = kwargs.get('lookup', None)
        self.precursor = kwargs.get('precursor', False)
        self.isobaric = kwargs.get('isobaric', False)
        self.rt_tol = kwargs.get('rttol', None)
        self.mz_tol = kwargs.get('mztol', None)
        self.mz_toltype = kwargs.get('mztoltype', None)

    def get_psms(self):
        """Creates iterator to write to new tsv. Contains input tsv
        lines plus quant data for these."""
        self.header, isob_header = prep.get_full_and_isobaric_headers(
            self.oldheader, self.quantdb, self.isobaric, self.precursor)
        self.psms = prep.generate_psms_quanted(self.quantdb, self.fn,
                                               isob_header, self.oldheader,
                                               self.isobaric, self.rt_tol,
                                               self.mz_tol, self.mz_toltype)
