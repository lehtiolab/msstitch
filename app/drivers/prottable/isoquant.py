from app.readers import tsv as reader
from app.actions.prottable import isoquant as preparation
from app.drivers.prottable.base import ProttableAddData


class AddIsobaricQuantDriver(ProttableAddData):
    outsuffix = '_isoq.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.quantfile = kwargs.get('quantfile', False)
        self.headertypes = ['isoquant']

    def initialize_input(self):
        super().initialize_input()
        self.quantheader = reader.get_tsv_header(self.quantfile)
        self.quantproteins = reader.generate_tsv_proteins(self.quantfile,
                                                          self.quantheader)

    def set_feature_generator(self):
        self.features = preparation.add_isoquant_data(self.in_proteins,
                                                      self.quantproteins,
                                                      self.quantpattern,
                                                      self.quantheader)
