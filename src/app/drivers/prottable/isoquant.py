from app.readers import tsv as reader
from app.actions.prottable import isoquant as preparation
from app.drivers.prottable.base import ProttableAddData, ProttableDriver
from app.drivers.options import prottable_options


class IsobaricQuantDriver(ProttableDriver):
    outsuffix = '_isoquant.txt'
    command = 'isoquant'
    commandhelp = ('Creates a protein quantification table from a PSM table '
                   'input with isobaric quant data.')

    def __init__(self):
        super().__init__()
        self.infiletype = 'TSV PSM table (MSGF+)'

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['proteincol',
                                                 'quantcolpattern']))

    def create_header(self):
        self.header = preparation.get_isoquant_header(self.quantfields)

    def initialize_input(self):
        psmheader = reader.get_tsv_header(self.fn)
        self.get_column_header_for_number(['proteincol'], psmheader)
        self.quantfields = reader.get_columns_by_pattern(psmheader,
                                                         self.quantcolpattern)
        self.in_psms = reader.generate_tsv_proteins(self.fn,
                                                    psmheader)

    def set_feature_generator(self):
        self.features = preparation.isobaric_quant_psms(self.in_psms,
                                                        self.proteincol,
                                                        self.quantfields)


class AddIsobaricQuantDriver(ProttableAddData):
    outsuffix = '_added_isoq.txt'
    command = 'addisoquant'
    commandhelp = ('Add isobaric quantification data from another '
                   'proteintable containing this.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['quantfile',
                                                 'quantacccolpattern',
                                                 'quantcolpattern'],
                                                prottable_options))

    def initialize_input(self):
        super().initialize_input()
        quantheader = reader.get_tsv_header(self.quantfile)
        self.quantfields = reader.get_cols_in_file(self.quantcolpattern,
                                                   quantheader)
        self.quantacc = reader.get_cols_in_file(self.quantacccolpattern,
                                                quantheader, single_col=True)
        self.quantfeatures = reader.generate_tsv_proteins(self.quantfile,
                                                          quantheader)

    def create_header(self):
        self.header = self.oldheader + self.quantfields

    def set_feature_generator(self):
        self.features = preparation.add_isoquant_data(self.in_proteins,
                                                      self.quantfeatures,
                                                      self.quantacc,
                                                      self.quantfields)
