from app.readers import tsv as reader
from app.writers import prottable as writers
from app.actions.headers import prottable as head
from app.drivers.pepprottable import PepProttableDriver


class ProttableDriver(PepProttableDriver):
    def __init__(self):
        super().__init__()
        self.infiletype = 'protein table'

    def create_header(self):
        self.headerfields = head.get_prottable_headerfields(self.headertypes,
                                                            self.lookup,
                                                            self.poolnames)
        self.header = head.generate_header(self.headerfields, self.oldheader,
                                           self.group_by_field)


class ProttableAddData(ProttableDriver):
    def initialize_input(self):
        self.oldheader = reader.get_tsv_header(self.fn)
        self.in_proteins = reader.generate_tsv_proteins(self.fn,
                                                        self.oldheader)


class ProttableMergeDriver(ProttableDriver):
    def __init__(self):
        super().__init__()
        self.group_by_field = True

    def initialize_input(self):
        self.headertypes = ['proteindata']
        for inflag, htype in zip([self.probability, self.fdr,
                                  self.pep, self.precursor,
                                  self.isobaric],
                                 ['probability', 'proteinfdr',
                                  'proteinpep', 'precursorquant', 'isoquant']):
            if inflag:
                self.headertypes.append(htype)
        self.poolnames = [x[0] for x in self.lookup.get_all_poolnames()]

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_prottable_with_na(self.header, self.features, outfn)
