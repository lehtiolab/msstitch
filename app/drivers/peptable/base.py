from app.drivers.pepprottable import PepProttableDriver
from app.actions.headers import peptable as head


class PeptableDriver(PepProttableDriver):
    def create_header(self):
        self.headerfields = head.get_peptable_headerfields(self.headertypes,
                                                           self.lookup,
                                                           self.poolnames)
        self.header = head.generate_header(self.headerfields,
                                           self.oldheader)


class PeptableMergeDriver(PeptableDriver):
    def initialize_input(self):
        self.headertypes = ['proteindata']
        for inflag, htype in zip([self.fdr, self.pep, self.precursorquant,
                                  self.isobaricquant, self.nopsms],
                                 ['peptidefdr', 'peptidepep', 'precursorquant',
                                  'isoquant', 'nopsms']):
            if inflag:
                self.headertypes.append(htype)
        self.poolnames = [x[0] for x in self.lookup.get_all_poolnames()]
