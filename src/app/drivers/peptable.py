import sys
from app.drivers.base import PepProttableDriver
from app.drivers.options import peptable_options

from app.readers import tsv as tsvreader

from app.dataformats import mzidtsv as psmh
from app.dataformats import peptable as peph 
from app.dataformats import prottable as proth

from app.actions.psmtable import isosummarize
from app.actions import psmtopeptable


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
            'keep_psms_na', 'minpeptidenr', 'totalprotfn'], peptable_options))

    def prepare(self):
        self.oldheader = tsvreader.get_tsv_header(self.fn)
        self.get_column_header_for_number(['spectracol'])
        self.scorecol = tsvreader.get_cols_in_file(self.scorecolpattern,
                                                   self.oldheader, True)
        self.precurquantcol = psmtopeptable.get_quantcols(self.precursorquantcolpattern,
                                                 self.oldheader, 'precur')


    def set_features(self):
        qpat = self.quantcolpattern if self.quantcolpattern else '[a-z]+[0-9]+plex_'
        header = [x for x in self.oldheader if x != psmh.HEADER_SPECFILE]
        try:
            isocols = tsvreader.get_columns_by_pattern(header, qpat)
        except RuntimeError:
            pass
        else:
            for col in isocols:
                header.pop(header.index(col))
        if self.precurquantcol:
            header = [peph.HEADER_AREA if x == self.precurquantcol
                      else x for x in header]
        header = [peph.HEADER_PEPTIDE, peph.HEADER_LINKED_PSMS] + [
                x for x in header if x != psmh.HEADER_PEPTIDE]
        switch_map = {old: new for old, new in zip(
            [psmh.HEADER_PEPTIDE, psmh.HEADER_PROTEIN, psmh.HEADER_PEPTIDE_Q],
            [peph.HEADER_PEPTIDE, peph.HEADER_PROTEINS, peph.HEADER_QVAL])}
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
            totalproteome, tpacc, tp_pepacc = False, False, False
            if self.totalprotfn:
                pep_tp_accs = [psmh.HEADER_MASTER_PROT, psmh.HEADER_SYMBOL,
                        psmh.HEADER_GENE, peph.HEADER_PROTEINS]
                totalphead = tsvreader.get_tsv_header(self.totalprotfn)
                totalproteome = tsvreader.generate_split_tsv_lines(self.totalprotfn, totalphead)
                totalpfield_found = False
                for tpacc, tp_pepacc in zip(proth.TPROT_HEADER_ACCS, pep_tp_accs):
                    if totalphead[0] == tpacc and tp_pepacc in self.header:
                        totalpfield_found = True
                        break
                if not totalpfield_found:
                    print('Could not find correct header field name in the total '
                            'proteome table passed. '
                            'Should be one of {}'.format(proth.TPROT_HEADER_ACCS))
                    sys.exit(1)
            nopsms = [isosummarize.get_no_psms_field(qf) for qf in quantcols]
            self.header = self.header + quantcols + nopsms + [proth.HEADER_NO_FULLQ_PSMS]
            peptides = isosummarize.get_isobaric_ratios(self.fn, self.oldheader, 
                    quantcols, denomcols, self.mediansweep, self.medianintensity,
                    self.median_or_avg, self.minint, peptides, self.header[0],
                    psmh.HEADER_PEPTIDE, totalproteome, tpacc, tp_pepacc,
                    self.logisoquant, self.mediannormalize, self.keepnapsms)
        if self.modelqvals:
            qix = self.header.index(peph.HEADER_QVAL) + 1
            self.header = self.header[:qix] + [peph.HEADER_QVAL_MODELED] + self.header[qix:]
            scorecol = tsvreader.get_cols_in_file(self.scorecolpattern,
                    self.oldheader, True)
            peptides = psmtopeptable.recalculate_qvals_linear_model(peptides, 
                    scorecol, self.qvalthreshold, self.minpeptidenr)
        self.features = peptides
