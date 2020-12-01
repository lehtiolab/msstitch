import sys

from app.drivers.base import PepProttableDriver

from app.drivers.options import prottable_options
from app.readers import tsv as tsvreader
from app.dataformats import prottable as prottabledata
from app.dataformats import peptable as peptabledata
from app.dataformats import mzidtsv as mzidtsvdata

from app.actions import proteins
from app.actions.psmtable import isosummarize


class ProttableDriver(PepProttableDriver):

    def set_options(self):
        super().set_options()
        options = self.define_options(['decoyfn', 'scorecolpattern', 'minlogscore',
            'quantcolpattern', 'minint', 'denomcols', 'denompatterns', 'mediansweep',
            'medianintensity', 'median_or_avg', 'logisoquant', 'mediannormalize',
            'keep_psms_na', 'precursor', 'psmfile'], prottable_options)
        self.options.update(options)

    def get_td_proteins_bestpep(self, theader, dheader):
        self.header = [self.headeraccfield] + prottabledata.PICKED_HEADER
        tscorecol = tsvreader.get_cols_in_file(self.scorecolpattern, theader, True)
        dscorecol = tsvreader.get_cols_in_file(self.scorecolpattern, dheader, True)
        tpeps = tsvreader.generate_split_tsv_lines(self.fn, theader)
        dpeps = tsvreader.generate_split_tsv_lines(self.decoyfn, dheader)
        targets = proteins.generate_bestpep_proteins(tpeps, tscorecol, 
                self.minlogscore, self.headeraccfield, self.fixedfeatcol)
        decoys = proteins.generate_bestpep_proteins(dpeps, dscorecol,
                self.minlogscore, self.headeraccfield, self.fixedfeatcol)
        return targets, decoys

    def get_quant(self, theader, features):
        if self.precursor:
            tpeps = tsvreader.generate_split_tsv_lines(self.fn, theader)
            self.header.append(prottabledata.HEADER_AREA)
            features = proteins.add_ms1_quant_from_top3_mzidtsv(features, 
                    tpeps, self.headeraccfield, self.fixedfeatcol)
        if self.quantcolpattern:
            psmheader = tsvreader.get_tsv_header(self.psmfile)
            denomcols = False
            if self.denomcols is not None:
                denomcols = [self.number_to_headerfield(col, psmheader)
                             for col in self.denomcols]
            elif self.denompatterns is not None:
                denomcolnrs = [tsvreader.get_columns_by_pattern(psmheader, pattern)
                               for pattern in self.denompatterns]
                denomcols = set([col for cols in denomcolnrs for col in cols])
            elif not self.mediansweep and not self.medianintensity:
                print('Must define either denominator column numbers '
                        'or regex pattterns to find them, or use median sweep, or '
                        'report median intensities.')
                sys.exit(1)
            elif self.medianintensity and self.mediannormalize:
                print('Cannot do median-centering on intensity values, exiting')
                sys.exit(1)
            quantcols = tsvreader.get_columns_by_pattern(psmheader, self.quantcolpattern)
            nopsms = [isosummarize.get_no_psms_field(qf) for qf in quantcols]
            self.header = self.header + quantcols + nopsms + [prottabledata.HEADER_NO_FULLQ_PSMS]
            features = isosummarize.get_isobaric_ratios(self.psmfile, psmheader,
                    quantcols, denomcols, self.mediansweep, self.medianintensity,
                    self.median_or_avg, self.minint, features, self.headeraccfield,
                    self.fixedfeatcol, False, False, False, self.logisoquant, self.mediannormalize,
                    self.keepnapsms)
        return features




class ProteinsDriver(ProttableDriver):
    command = 'proteins'
    commandhelp = 'Create a protein table from peptides'
    outsuffix = '_proteins.tsv'
    headeraccfield = prottabledata.HEADER_PROTEIN
    fixedfeatcol = peptabledata.HEADER_MASTERPROTEINS

    def set_features(self):
        theader = tsvreader.get_tsv_header(self.fn)
        dheader = tsvreader.get_tsv_header(self.decoyfn)
        targets, decoys = self.get_td_proteins_bestpep(theader, dheader)
        features = proteins.generate_protein_fdr(targets, decoys, self.headeraccfield)
        self.features = self.get_quant(theader, features)


class GenesDriver(ProttableDriver):
    command = 'genes'
    commandhelp = 'Create a gene table from peptides'
    outsuffix = '_genes.tsv'
    headeraccfield = prottabledata.HEADER_GENENAME
    fixedfeatcol = mzidtsvdata.HEADER_SYMBOL

    def set_options(self):
        super().set_options()
        options = self.define_options(['fastadelim', 'genefield', 't_fasta', 
            'd_fasta'], prottable_options)
        self.options.update(options)

    def set_features(self):
        theader = tsvreader.get_tsv_header(self.fn)
        dheader = tsvreader.get_tsv_header(self.decoyfn)
        targets, decoys = self.get_td_proteins_bestpep(theader, dheader)
        fastadelim, genefield = self.get_fastadelim_genefield(self.fastadelim,
                                                             self.genefield)
        features = proteins.generate_pick_fdr(
           targets, decoys, self.t_fasta, self.d_fasta, 'fasta', self.headeraccfield,
           fastadelim, genefield)
        self.features = self.get_quant(theader, features)


class ENSGDriver(GenesDriver):
    command = 'ensg'
    commandhelp = 'Create an ENSG table from peptides'
    outsuffix = '_ensg.tsv'
    headeraccfield = prottabledata.HEADER_GENEID
    fixedfeatcol = mzidtsvdata.HEADER_GENE
