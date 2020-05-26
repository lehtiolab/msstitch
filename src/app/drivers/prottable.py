from app.drivers.base import PepProttableDriver

from app.drivers.options import prottable_options
from app.readers import tsv as tsvreader
from app.dataformats import prottable as prottabledata
from app.dataformats import peptable as peptabledata
from app.dataformats import mzidtsv as mzidtsvdata

from app.actions import proteins
from app.actions.psmtable import isosummarize


class ProttableDriver(PepProttableDriver):
    mediannormalize = False # FIXME remove when done

    def set_options(self):
        super().set_options()
        options = self.define_options(['decoyfn', 'scorecolpattern', 'minlogscore',
            'quantcolpattern', 'minint', 'denomcols', 'denompatterns', 
            'precursor', 'psmfile'], prottable_options)
        self.options.update(options)

    def get_td_proteins_bestpep(self, theader, dheader):
        self.header = [self.headeraccfield] + prottabledata.PICKED_HEADER
        tscorecol = tsvreader.get_cols_in_file(self.scorecolpattern, theader, True)
        dscorecol = tsvreader.get_cols_in_file(self.scorecolpattern, dheader, True)
        tpeps = tsvreader.generate_tsv_psms(self.fn, theader)
        dpeps = tsvreader.generate_tsv_psms(self.decoyfn, dheader)
        targets = proteins.generate_bestpep_proteins(tpeps, tscorecol, 
                self.minlogscore, self.headeraccfield, self.featcol)
        decoys = proteins.generate_bestpep_proteins(dpeps, dscorecol,
                self.minlogscore, self.headeraccfield, self.featcol)
        return targets, decoys

    def get_quant(self, theader, features):
        if self.precursor:
            tpeps = tsvreader.generate_tsv_psms(self.fn, theader)
            self.header.append(prottabledata.HEADER_AREA)
            features = proteins.add_ms1_quant_from_top3_mzidtsv(features, 
                    tpeps, self.headeraccfield, self.featcol)
        if self.quantcolpattern:
            psmheader = tsvreader.get_tsv_header(self.psmfile)
            if self.denomcols is not None:
                denomcols = [self.number_to_headerfield(col, psmheader)
                             for col in self.denomcols]
            elif self.denompatterns is not None:
                denomcolnrs = [tsvreader.get_columns_by_pattern(psmheader, pattern)
                               for pattern in self.denompatterns]
                denomcols = set([col for cols in denomcolnrs for col in cols])
            else:
                raise RuntimeError('Must define either denominator column numbers '
                                   'or regex pattterns to find them')
            quantcols = tsvreader.get_columns_by_pattern(psmheader, self.quantcolpattern)
            nopsms = [isosummarize.get_no_psms_field(qf) for qf in quantcols]
            self.header = self.header + quantcols + nopsms
            features = isosummarize.get_isobaric_ratios(self.psmfile, psmheader,
                                                 quantcols, denomcols, self.minint,
                                                 features, self.headeraccfield,
                                                 self.featcol, self.mediannormalize)
        return features




class ProteinsDriver(ProttableDriver):
    command = 'proteins'
    commandhelp = 'Create a protein table from peptides'
    outsuffix = '_proteins.tsv'
    headeraccfield = prottabledata.HEADER_PROTEIN
    featcol = peptabledata.HEADER_MASTERPROTEINS

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
    featcol = mzidtsvdata.HEADER_SYMBOL

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
    featcol = mzidtsvdata.HEADER_GENE


# TODO create a result driver? For what?
