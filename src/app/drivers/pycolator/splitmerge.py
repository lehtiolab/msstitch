from app.drivers.pycolator import base
from app.actions.pycolator import splitmerge as preparation
from app.readers import pycolator as readers
from app.drivers.options import pycolator_options


class SplitDriver(base.PycolatorDriver):
    outfile = None

    def run(self):
        self.set_filter_types()
        for filter_type, suffix in self.filter_types:
            self.prepare()
            self.set_features(filter_type)
            self.outsuffix = suffix
            self.write()
        self.finish()

    def set_options(self):
        """Since splitdriver splits into multiple files we cannot set an
        output file"""
        super().set_options()
        del(self.options['-o'])

    def set_features(self, filter_type):
        """Calls splitter to split percolator output into target/decoy
        elements.
        Writes two new xml files with features. Currently only psms and
        peptides. Proteins not here, since one cannot do protein inference
        before having merged and remapped multifraction data anyway.
        """
        elements_to_split = {'psm': self.allpsms, 'peptide': self.allpeps}
        self.features = self.splitfunc(elements_to_split, self.ns, filter_type)


class SplitProteinDriver(SplitDriver):
    """Using --protheaders "novel_ENSP" "ENSP|variant_"  will result in two output files,
    the first containing all PSM/peptides containing at least one mapping with a header novel_ENSP,
    and the second containing all PSMs/peptides containing at least one mapping with header variant_,
    but no PSMs/peptides which map to ENSP headers.
    """
    command = 'splitprotein'
    commandhelp = ('Splits input XML into multiple files depending based on '
                   'the protein headers specified. Each header class gets '
                   'its own output file')

    def set_filter_types(self):
        maxdigits = len(str(len(self.protheaders)))
        self.filter_types = [(headers, '_h{i:0{dig}d}.xml'.format(
            i=ix, dig=maxdigits))
            for ix, headers in enumerate(self.protheaders)]

    def set_features(self, filter_type):
        self.splitfunc = preparation.split_protein_header_id_type
        super().set_features(filter_type)

    def set_options(self):
        super().set_options()
        options = self.define_options(['protheaders'], pycolator_options)
        self.options.update(options)
