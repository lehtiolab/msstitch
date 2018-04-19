from app.drivers.pepprottable import PepProttableDriver
from app.writers import prottable as writers
from app.actions.headers import peptable as head
from app.readers import tsv as reader


class PeptableDriver(PepProttableDriver):
    def __init__(self):
        super().__init__()
        self.infiletype = 'peptide table'

    def create_header(self):
        self.headerfields = head.get_peptable_headerfields(self.headertypes,
                                                           self.lookup,
                                                           self.poolnames)
        self.header = head.generate_header(self.headerfields, self.oldheader,
                                           self.group_by_field)


class PeptableAddData(PeptableDriver):
    def initialize_input(self):
        self.oldheader = reader.get_tsv_header(self.fn)
        self.in_peptides = reader.generate_tsv_proteins(self.fn,
                                                        self.oldheader)


class PeptableMergeDriver(PeptableDriver):
    def __init__(self):
        super().__init__()
        self.group_by_field = True

    def initialize_input(self):
        self.headertypes = ['proteindata']
        for inflag, htype in zip([self.fdr, self.pep, self.precursor,
                                  self.isobaric],
                                 ['peptidefdr', 'peptidepep', 'precursorquant',
                                  'isoquant']):
            if inflag:
                self.headertypes.append(htype)
        self.poolnames = [x[0] for x in self.lookup.get_all_poolnames()]

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_prottable_with_na(self.header, self.features, outfn)
