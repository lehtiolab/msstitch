from app.drivers.base import BaseDriver
from app.writers import prottable as writers
from app.readers import tsv as reader
from app.preparation import prottable as preparation
from app.lookups.sqlite.proteingroups import ProteinGroupProteinTableDB


class AddProteinInfoDriver(BaseDriver):
    outsuffix = '_proteindata.txt'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.fasta = kwargs.get('fasta')
        self.lookup = kwargs.get('protgroupdb')

    def run(self):
        self.oldheader = reader.get_tsv_header(self.fn)
        self.set_protein_generator()
        self.write()
        self.finish()

    def write(self):
        outfn = self.create_outfilepath(self.fn, self.outsuffix)
        writers.write_prottable(self.header, self.proteins, outfn)

    def set_protein_generator(self):
        protgroupdb = ProteinGroupProteinTableDB(self.lookup)
        preparation.collect_descriptions(protgroupdb, self.fasta)
        self.header = preparation.get_header_with_proteindata(self.oldheader)
        proteins = reader.generate_tsv_proteins(self.fn, self.oldheader)
        self.proteins = preparation.add_protein_data(proteins,
                                                     self.oldheader,
                                                     protgroupdb)
