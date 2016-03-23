from app.readers import tsv as reader
from app.actions.peptable import isoquant as preparation
from app.drivers.peptable.base import PeptableAddData
from app.drivers.options import peptable_options


class AddIsobaricQuantDriver(PeptableAddData):
    outsuffix = '_isoq.txt'
    command = 'isoquant'
    commandhelp = ('Add isobaric quantification data from another '
                   'peptidetable containing this.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['quantfile',
                                                 'quantcolpattern',
                                                 'quantacccolpattern'],
                                                peptable_options))
        self.options['--isobquantcolpattern'].update({'required': True})

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.headertypes = ['isoquant']

    def initialize_input(self):
        super().initialize_input()
        quantheader = reader.get_tsv_header(self.quantfile)
        self.quantfields = reader.get_cols_in_file(self.quantcolpattern,
                                                   quantheader)
        self.quantacc = reader.get_cols_in_file(self.quantacccolpattern,
                                                quantheader, single_col=True)
        self.quantpeptides = reader.generate_tsv_proteins(self.quantfile,
                                                          quantheader)

    def create_header(self):
        self.header = self.oldheader + self.quantfields

    def set_feature_generator(self):
        self.features = preparation.add_isoquant_data(self.in_peptides,
                                                      self.quantpeptides,
                                                      self.quantacc,
                                                      self.quantfields)
