from app.drivers.base import PepProttableDriver
from app.drivers.options import peptable_options

from app.readers import tsv as tsvreader

from app.dataformats import mzidtsv as mzidtsvdata
from app.dataformats import peptable as peptabledata 
from app.dataformats.prottable import HEADER_NO_FULLQ_PSMS

from app.actions.psmtable import isosummarize
from app.actions import psmtopeptable


# TODO:
# PSM table to pep/prot quant without ratios, so "median intensity" only
# But, also retain a peptide table without quant, so DEqMS can fill in normalized


class CreatePeptableDriver(PepProttableDriver):
    """Creates unique peptide table from MzidTSV table. Will not change
    q-values/PEP after filtration, so these should have been calculated
    for the peptides beforehand or rescored afterwards.
    """
    outsuffix = '_peptable.tsv'
    command = 'peptides'
    commandhelp = ('Create peptide table from PSM TSV input, uses best '
                   'scoring PSM for each peptide and strips PSM '
                   'isobaric quant information. Retains MS1 quant info if '
                   'specified, by taking the highest precursor quant value '
                   'for a peptide. Optionally add isobaric quant and  modeled '
                   'q-values')

    def __init__(self):
        super().__init__()
        self.infiletype = 'TSV PSM table (MSGF+)'

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['spectracol', 'scorecolpattern', 
            'quantcolpattern', 'precursorquantcolpattern', 'minint', 'denomcols',
            'denompatterns', 'mediansweep', 'medianintensity', 'median_or_avg',
            'logisoquant', 'mediannormalize', 'modelqvals', 'qvalthreshold',
            'minpeptidenr'], peptable_options))

    def prepare(self):
        self.oldheader = tsvreader.get_tsv_header(self.fn)
        self.get_column_header_for_number(['spectracol'])
        self.scorecol = tsvreader.get_cols_in_file(self.scorecolpattern,
                                                   self.oldheader, True)
        self.precurquantcol = psmtopeptable.get_quantcols(self.precursorquantcolpattern,
                                                 self.oldheader, 'precur')


    def set_features(self):
        qpat = self.quantcolpattern if self.quantcolpattern else '[a-z]+[0-9]+plex_'
        header = [x for x in self.oldheader if x != mzidtsvdata.HEADER_SPECFILE]
        try:
            isocols = tsvreader.get_columns_by_pattern(header, qpat)
        except RuntimeError:
            pass
        else:
            for col in isocols:
                header.pop(header.index(col))
        if self.precurquantcol:
            header = [peptabledata.HEADER_AREA if x == self.precurquantcol
                      else x for x in header]
        header = [peptabledata.HEADER_PEPTIDE, peptabledata.HEADER_LINKED_PSMS] + [
                x for x in header if x != mzidtsvdata.HEADER_PEPTIDE]
        switch_map = {old: new for old, new in zip(
            [mzidtsvdata.HEADER_PEPTIDE, mzidtsvdata.HEADER_PROTEIN, mzidtsvdata.HEADER_PEPTIDE_Q],
            [peptabledata.HEADER_PEPTIDE, peptabledata.HEADER_PROTEINS, peptabledata.HEADER_QVAL])}
        self.header = [switch_map[field] if field in switch_map else field
                for field in header]
        peptides = psmtopeptable.generate_peptides(self.fn, self.oldheader,
                switch_map, self.scorecol, self.precurquantcol, self.spectracol)
        # Remove quant data if not specified any way to summarize
        if self.quantcolpattern and any([self.denomcols, self.denompatterns,
                self.mediansweep, self.medianintensity]):
            denomcols = False
            if self.denomcols is not None:
                denomcols = [self.number_to_headerfield(col, self.oldheader)
                             for col in self.denomcols]
            elif self.denompatterns is not None:
                denomcolnrs = [tsvreader.get_columns_by_pattern(self.oldheader, pattern)
                               for pattern in self.denompatterns]
                denomcols = set([col for cols in denomcolnrs for col in cols])
            quantcols = tsvreader.get_columns_by_pattern(self.oldheader,
                                                   self.quantcolpattern)
            nopsms = [isosummarize.get_no_psms_field(qf) for qf in quantcols]
            self.header = self.header + quantcols + nopsms + [HEADER_NO_FULLQ_PSMS]
            peptides = isosummarize.get_isobaric_ratios(self.fn, self.oldheader, 
                    quantcols, denomcols, self.mediansweep, self.medianintensity,
                    self.median_or_avg, self.minint, peptides, self.header[0],
                    mzidtsvdata.HEADER_PEPTIDE, self.logisoquant, self.mediannormalize)
        if self.modelqvals:
            qix = self.header.index(peptabledata.HEADER_QVAL) + 1
            self.header = self.header[:qix] + [peptabledata.HEADER_QVAL_MODELED] + self.header[qix:]
            scorecol = tsvreader.get_cols_in_file(self.scorecolpattern,
                    self.oldheader, True)
            peptides = psmtopeptable.recalculate_qvals_linear_model(peptides, 
                    scorecol, self.qvalthreshold, self.minpeptidenr)
        self.features = peptides
