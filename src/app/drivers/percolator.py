from app.drivers import base
from app.drivers.options import percolator_options

from app.actions.percolator import split
from app.actions.percolator import filters


class FilterSequences(base.PercolatorDriver):
    """This class processes multiple percolator runs from fractions and
    filters out first peptides that are found in a specified searchspace. Then
    it keeps the remaining best scoring unique peptides."""
    outsuffix = '_filtered.xml'
    command = 'filterperco'
    lookuptype = 'searchspace'
    commandhelp = ('Filters out peptides that match to a protein in a FASTA '
                   'file. Needs a lookup to match peptides to quickly.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['fullprotein', 'deamidate', 
            'fasta', 'minlength', 'lookupfn', 'forcetryp', 'falloff'],
            percolator_options))

    def set_features(self):
        if self.fullprotein:
            self.features = {
                'peptide': filters.filter_whole_proteins(self.allpeps,
                                                             self.fasta,
                                                             self.lookup,
                                                             'pep',
                                                             self.ns,
                                                             self.deamidate,
                                                             self.minlength,
                                                             self.forcetryp),
                'psm': filters.filter_whole_proteins(self.allpsms,
                                                         self.fasta,
                                                         self.lookup,
                                                         'psm',
                                                         self.ns,
                                                         self.deamidate,
                                                         self.minlength,
                                                         self.forcetryp)
            }
        else:
            self.features = {
                'peptide': filters.filter_known_searchspace(self.allpeps,
                                                                'pep',
                                                                self.lookup,
                                                                self.ns,
                                                                self.falloff,
                                                                self.deamidate),
                'psm': filters.filter_known_searchspace(self.allpsms,
                                                            'psm',
                                                            self.lookup,
                                                            self.ns,
                                                            self.falloff,
                                                            self.deamidate),
            }


class SplitDriver(base.PercolatorDriver):
    outfile = None

    def run(self):
        self.set_filter_types()
        for filter_type, suffix in self.filter_types:
            self.prepare()
            self.set_features(filter_type)
            self.outsuffix = suffix
            self.write()

    def set_options(self):
        """Since splitdriver splits into multiple files we cannot set an
        output file"""
        super().set_options()
        del(self.options['outfile'])

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
    command = 'splitperco'
    commandhelp = ('Splits input XML into multiple files depending based on '
                   'the protein headers specified. Each header class gets '
                   'its own output file')

    def set_filter_types(self):
        maxdigits = len(str(len(self.protheaders)))
        self.filter_types = [(headers, '_h{i:0{dig}d}.xml'.format(
            i=ix, dig=maxdigits))
            for ix, headers in enumerate(self.protheaders)]

    def set_features(self, filter_type):
        self.splitfunc = split.split_protein_header_id_type
        super().set_features(filter_type)

    def set_options(self):
        super().set_options()
        options = self.define_options(['protheaders'], percolator_options)
        self.options.update(options)
