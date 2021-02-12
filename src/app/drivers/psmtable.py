import os
import sys
from hashlib import md5
from itertools import chain

from app.drivers.options import psmtable_options
from app.drivers.base import PSMDriver

from app.readers import tsv as tsvreader
from app.readers import mzidplus as mzidreader
from app.dataformats import mzidtsv as psmhead
from app.dataformats.prottable import HEADER_NO_FULLQ_PSMS

from app.actions.psmtable import splitmerge as splitmerge
from app.actions.psmtable import refine as refine
from app.actions.psmtable import percolator as perco
from app.actions.psmtable import filtering
from app.actions.psmtable import isosummarize

from app.dataformats import prottable as prottabledata

from app.writers import tsv as writer


class TSVConcatenateDriver(PSMDriver):
    """Concatenates TSVs, takes care of headers"""
    outsuffix = '_concat.tsv'
    command = 'concat'
    commandhelp = ('Merges multiple TSV tables of MSGF+ output.'
                   'Make sure headers are same in all files.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['multifiles'],
                                                psmtable_options))

    def prepare(self):
        # do not read PSMs, multiple files passed and they will be checked
        # if headers matched
        self.first_infile = self.fn[0]
        self.oldheader = tsvreader.get_tsv_header(self.first_infile)

    def set_features(self):
        self.header = self.oldheader
        self.psms = splitmerge.merge_mzidtsvs(self.fn, self.header)


class TSVSplitDriver(PSMDriver):
    """Splits MSGF PSM table on contents of certain column. Each
    row in file is piped to an output file. Which output files
    the row is written to depends on the contents of the selected
    column"""
    outsuffix = '_split.tsv'
    command = 'split'
    commandhelp = 'Splits an MSGF TSV PSM table into multiple new tables'

    def set_options(self):
        super().set_options()
        options = self.define_options(['splitcol'], psmtable_options)
        self.options.update(options)

    def set_features(self):
        self.header = self.oldheader[:]
        splitfield = splitmerge.get_splitfield(self.oldheader, self.splitcol)
        self.psms = splitmerge.generate_psms_split(self.oldpsms, splitfield)

    def write(self):
        base_outfile = os.path.join(self.outdir, '{}.tsv')
        writer.write_multi_mzidtsv(self.header, self.oldheader, self.psms,
                                    base_outfile)


class Perco2PSMDriver(PSMDriver):
    """
    Adds percolator data from mzid file to table.
    """
    outsuffix = '_fdr.tsv'
    command = 'perco2psm'
    commandhelp = (
            'Calculates FDR from percolator output and adds FDR and percolator '
            ' to the corresponding TSV PSM tables. FDR calculation method is TD-'
            'competition.'
            )

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['multifiles', 'mzidfns', 'percofn',
            'filtpep', 'filtpsm'], psmtable_options))

    def prepare(self):
        # multiple PSM tables passed so do not read here, match with mzid
        pass

    def set_features(self):
        self.percopsms = perco.calculate_target_decoy_competition(self.percofn)
                
    def write(self):
        for psmfn, mzidfn in zip(self.fn, self.mzidfns):
            oldheader = tsvreader.get_tsv_header(psmfn)
            header = perco.get_header_with_percolator(oldheader)
            outfn = self.create_outfilepath(psmfn, self.outsuffix)
            mzns = mzidreader.get_mzid_namespace(mzidfn)
            mzidsr = mzidreader.mzid_spec_result_generator(mzidfn, mzns)
            psms = tsvreader.generate_split_tsv_lines(psmfn, oldheader)
            psms_perco = perco.add_fdr_to_mzidtsv(psms, mzidsr, mzns,
                    self.percopsms)
            if self.filtpsm:
                psms_perco = filtering.filter_psms_conf(psms_perco, psmhead.HEADER_PSMQ,
                        self.filtpsm, True)
            if self.filtpep:
                psms_perco = filtering.filter_psms_conf(psms_perco, psmhead.HEADER_PEPTIDE_Q,
                        self.filtpep, True)
            writer.write_tsv(header, psms_perco, outfn)


class ConfidenceFilterDriver(PSMDriver):
    outsuffix = '_filtconf.txt'
    command = 'conffilt'
    commandhelp = 'Filters PSMs by their confidence level. '

    def set_options(self):
        super().set_options()
        options = self.define_options(['confcol', 'confpattern', 'conflvl',
                                       'conftype', 'unroll'], psmtable_options)
        self.options.update(options)

    def parse_input(self, **kwargs):
        super().parse_input(**kwargs)
        self.lowerbetter = self.conftype == 'lower'

    def set_features(self):
        self.header = self.oldheader[:]
        if self.confpattern:
            confkey = tsvreader.get_cols_in_file(self.confpattern,
                                                 self.header, True)
        elif self.confcol:
            confkey = self.header[int(self.confcol) - 1]
        else:
            print('Must define either --confcol or --confcolpattern')
            sys.exit(1)
        self.psms = filtering.filter_psms_conf(self.oldpsms,
                                       confkey,
                                       self.conflvl,
                                       self.lowerbetter)


class PSMTableRefineDriver(PSMDriver):
    # gene, quant, pg, spectra
    outsuffix = '_refined.txt'
    lookuptype = 'psm'
    command = 'psmtable'
    commandhelp = ('Add columns to mzidtsv with genes, quantification, protein groups, '
                   'etc. which are stored in a lookup specified with --dbfile')

    def set_options(self):
        super().set_options()
        options = self.define_options(['oldpsmfile', 'lookupfn', 'precursor', 'isobaric',
            'minpurity', 'unroll', 'spectracol', 'addbioset', 'addmiscleav', 'genes',
            'proteingroup', 'fasta', 'genefield', 'fastadelim'], psmtable_options)
        self.options.update(options)

    def set_features(self):
        """Creates iterator to write to new tsv. Contains input tsv
        lines plus quant data for these."""
        # First prepare the data, read PSM table to SQLite
        specfncolnr = int(self.spectracol) - 1
        specfncol = self.oldheader[specfncolnr]
        fastadelim, genefield = self.get_fastadelim_genefield(self.fastadelim,
                                                              self.genefield)

        if self.fasta:
            fasta_md5 = refine.get_fasta_md5(self.fasta)
        else:
            fasta_md5 = False

        # If appending to previously refined PSM table, reuse DB and shift rows
        if self.oldpsmfile:
            oldfasta_md5 = self.lookup.get_fasta_md5()
            if fasta_md5 != oldfasta_md5:
                print('WARNING, FASTA database used in old PSM table differs '
                        'from the passed database (or this cannot be determined '
                        'due to version differences), this may cause problems, as '
                        'msstitch will use the old database for PSM annotation.')
            shiftrows = self.lookup.get_highest_rownr() + 1
            proteins = set([x for x in self.lookup.get_protids()])
            self.lookup.drop_psm_indices()
        else:
            shiftrows = 0

        if self.proteingroup:
            if not fasta_md5 and not oldfasta_md5:
                # In case of old Fasta already stored it will be fine to protein group
                print('Cannot create protein group without supplying FASTA search '
                        'database file')
                sys.exit(1)
            self.tabletypes.append('proteingroup')
            self.lookup.drop_pgroup_tables()
        self.lookup.add_tables(self.tabletypes)

        # Need to place this here since we cannot store before having done add tables, but that
        # has to be done after getting proteingroup knowledge, which depends on knowledge of 
        # having passed an oldpsmfile (because of oldfasta_md5):
        if not self.oldpsmfile:
            proteins = refine.store_proteins_descriptions(self.lookup, self.fasta,
                    fasta_md5, self.fn, self.oldheader, fastadelim, genefield)

        refine.create_psm_lookup(self.fn, self.oldheader, proteins, self.lookup, 
                shiftrows, self.unroll, specfncol, fastadelim, genefield)
        isob_header = [x[0] for x in self.lookup.get_all_quantmaps()] if self.isobaric else False
        self.header = refine.create_header(self.oldheader, self.genes, 
                self.proteingroup, self.precursor, isob_header, self.addbioset, 
                self.addmiscleav, specfncolnr)
        psms = self.oldpsms
        # Now pass PSMs through multiple generators to add info
        if self.genes:
            psms = refine.add_genes_to_psm_table(psms, self.lookup)
        if self.isobaric or self.precursor:
            psms = refine.generate_psms_quanted(self.lookup, shiftrows, psms,
                    isob_header, self.isobaric, self.precursor, self.min_purity)
        psms = refine.generate_psms_spectradata(self.lookup, shiftrows, 
                psms, self.addbioset, self.addmiscleav)
        if self.oldpsmfile:
            prevheader = tsvreader.get_tsv_header(self.oldpsmfile)
            previouspsms = tsvreader.generate_split_tsv_lines(self.oldpsmfile, prevheader)
            psms = chain(previouspsms, psms)
        # Enforce proteingroup last, since it has to come AFTER the chaining of old + new PSMs
        # In theory you could do it before, but that makes no sense since in a big experiment you
        # also do not map PSMs to genes differently? If that is needed, you have to run multiple
        # experiments.
        if self.proteingroup:
            refine.build_proteingroup_db(self.lookup)
            psms = refine.generate_psms_with_proteingroups(psms, self.lookup, specfncol, self.unroll)
        self.psms = psms
        

class IsoSummarizeDriver(PSMDriver):
    outsuffix = '_ratio_isobaric.txt'
    command = 'isosummarize'
    commandhelp = ('Produce isobaric summarized data from a '
                   'PSM table containing raw intensities.')

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['quantcolpattern', 'denompatterns',
            'denomcols', 'mediansweep', 'medianintensity', 'median_or_avg',
            'keep_psms_na', 'minint', 'featcol', 'logisoquant', 
            'mediannormalize'], psmtable_options))

    def set_features(self):
        denomcols = False
        if self.denomcols is not None:
            denomcols = [self.number_to_headerfield(col, self.oldheader)
                         for col in self.denomcols]
        elif self.denompatterns is not None:
            denomcolnrs = [tsvreader.get_columns_by_pattern(self.oldheader, pattern)
                           for pattern in self.denompatterns]
            denomcols = set([col for cols in denomcolnrs for col in cols])
        elif not self.mediansweep and not self.medianintensity:
            raise RuntimeError('Must define either denominator column numbers '
                               'or regex pattterns to find them')
        quantcols = tsvreader.get_columns_by_pattern(self.oldheader,
                                               self.quantcolpattern)
        nopsms = [isosummarize.get_no_psms_field(qf) for qf in quantcols]
        if self.featcol:
            self.get_column_header_for_number(['featcol'], self.oldheader)
            self.header = [self.featcol] + quantcols + nopsms + [HEADER_NO_FULLQ_PSMS]
        else:
            self.header = (self.oldheader +
                           ['ratio_{}'.format(x) for x in quantcols])
        self.psms = isosummarize.get_isobaric_ratios(self.fn, self.oldheader,
                quantcols, denomcols, self.mediansweep, self.medianintensity,
                self.median_or_avg, self.minint, False, False, self.featcol,
                False, False, False, self.logisoquant, self.mediannormalize, self.keepnapsms)


class DeleteSetDriver(PSMDriver):
    outsuffix = '_deletedset.txt'
    command = 'deletesets'
    lookuptype = 'psm'
    commandhelp = """Remove sample sets from an existing PSM table and its lookup file"""

    def set_options(self):
        super().set_options()
        self.options.update(self.define_options(['lookupfn', 'setnames'], psmtable_options))
        self.options['lookupfn'].update({'required': False, 'default': None})

    def set_features(self):
        # In rare cases, a lookup is not used and we just filter the PSM table. I.e when
        # making a new lookup from a new+old PSM table.
        if self.lookup:
            self.lookup.delete_sample_set_shift_rows(self.setnames)
        self.header = self.oldheader
        self.psms = filtering.filter_psms_remove_set(self.oldpsms, self.setnames)
