from app.readers import tsv as reader
from app.actions.peptable import isoquant as preparation
from app.drivers.peptable.base import PeptableAddData


class AddIsobaricQuantDriver(PeptableAddData):
    outsuffix = '_isoq.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.quantfile = kwargs.get('quantfile', False)
        self.quantpattern = kwargs.get('quantcolpattern', False)
        self.quantacccol = kwargs.get('quantacccolpattern', False)
        self.headertypes = ['isoquant']

    def initialize_input(self):
        super().initialize_input()
        quantheader = reader.get_tsv_header(self.quantfile)
        self.quantfields = reader.get_cols_in_file(self.quantpattern, quantheader)
        self.quantacc = reader.get_cols_in_file(self.quantacccol, quantheader,
                                                single_col=True)
        self.quantpeptides = reader.generate_tsv_proteins(self.quantfile,
                                                          quantheader)

    def create_header(self):
        self.header = self.oldheader + self.quantfields

    def set_feature_generator(self):
        self.features = preparation.add_isoquant_data(self.in_peptides,
                                                      self.quantpeptides,
                                                      self.quantacc,
                                                      self.quantfields)
