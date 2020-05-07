import os
import re
from lxml import etree
from statistics import median

from app.dataformats import mzidtsv as constants
from tests.integration import basetests


class TestAddPSMData(basetests.MzidTSVBaseTest):
    command = 'specdata'
    suffix = '_spectradata.tsv'
    infilename = 'target.tsv'

    def test_addspec_miscleav_bioset(self):
        options = ['--dbfile', self.dbfile, '--spectracol', '1', '--addmiscleav', '--addbioset']
        self.run_command(options)
        sql = ('SELECT pr.rownr, bs.set_name, sp.retention_time, '
               'iit.ion_injection_time, im.ion_mobility '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'JOIN mzml AS sp USING(spectra_id) '
               'JOIN ioninjtime AS iit USING(spectra_id) '
               'JOIN ionmob AS im USING(spectra_id) '
               'JOIN mzmlfiles USING(mzmlfile_id) '
               'JOIN biosets AS bs USING(set_id) '
               'ORDER BY pr.rownr')
        fields = ['Biological set', 'Retention time(min)',
                  'Ion injection time(ms)', 'Ion mobility(Vs/cm2)']
        expected_values = self.process_dbvalues_both(self.dbfile, sql, [],
                                                     [1, 2, 3, 4], fields)
        self.check_results_sql(fields, self.rowify(expected_values))
        for val, exp in zip(self.get_values(['missed_cleavage']), self.get_values(['Peptide'], self.infile[0])):
            exp = re.sub('[0-9\+\.]', '', exp[0][1])[:-1]
            self.assertEqual(int(val[0][1]), exp.count('K') + exp.count('R') - exp.count('KP') - exp.count('RP'))

    def test_ionmobility(self):
        self.infilename = 'few_spec_timstof.tsv'
        self.dbfile = os.path.join(self.fixdir, 'td_psms_timstof.sqlite')
        options = ['--dbfile', self.dbfile, '--spectracol', '1', '--addmiscleav', '--addbioset']
        self.run_command(options)
        sql = ('SELECT pr.rownr, bs.set_name, sp.retention_time, '
               'iit.ion_injection_time, im.ion_mobility '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'JOIN mzml AS sp USING(spectra_id) '
               'JOIN ioninjtime AS iit USING(spectra_id) '
               'JOIN ionmob AS im USING(spectra_id) '
               'JOIN mzmlfiles USING(mzmlfile_id) '
               'JOIN biosets AS bs USING(set_id) '
               'ORDER BY pr.rownr')
        fields = ['Biological set', 'Retention time(min)',
                  'Ion injection time(ms)', 'Ion mobility(Vs/cm2)']
        expected_values = self.process_dbvalues_both(self.dbfile, sql, [],
                                                     [1, 2, 3, 4], fields)
        self.check_results_sql(fields, self.rowify(expected_values))
        for val, exp in zip(self.get_values(['missed_cleavage']), self.get_values(['Peptide'], self.infile[0])):
            exp = re.sub('[0-9\+\.]', '', exp[0][1])[:-1]
            self.assertEqual(int(val[0][1]), exp.count('K') + exp.count('R') - exp.count('KP') - exp.count('RP'))


class TestQuantTSV(basetests.MzidTSVBaseTest):
    command = 'quant'
    suffix = '_quant.tsv'
    infilename = 'target.tsv'

    def test_quanttsv_isobaric(self):
        options = ['--dbfile', self.dbfile, '--isobaric']
        self.run_command(options)
        sql = ('SELECT pr.rownr, ic.channel_name, iq.intensity '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'JOIN isobaric_quant AS iq USING(spectra_id) '
               'JOIN isobaric_channels AS ic USING(channel_id)')
        expected_values = self.get_values_from_db(self.dbfile, sql)
        fields = ['tmt10plex_{}'.format(ch) for ch in ['126', '127N', '127C',
                                                       '128N', '128C', '129N',
                                                       '129C', '130N', '130C',
                                                       '131']]
        self.check_results_sql(fields, self.rowify(expected_values))

    def test_quanttsv_precursor(self):
        options = ['--dbfile', self.dbfile, '--precursor']
        self.run_command(options)
        sql = ('SELECT pr.rownr, pq.intensity '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'LEFT OUTER JOIN ms1_align USING(spectra_id) '
               'LEFT OUTER JOIN ms1_quant AS pq USING(feature_id)')
        expected_values = self.get_values_from_db(self.dbfile, sql)
        self.check_results_sql(['MS1 area'], self.rowify(expected_values))

    def test_quanttsv_both(self):
        options = ['--dbfile', self.dbfile, '--isobaric', '--precursor']
        self.run_command(options)
        sql = ('SELECT pr.rownr, ic.channel_name, iq.intensity, pq.intensity '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'JOIN isobaric_quant AS iq USING(spectra_id) '
               'JOIN isobaric_channels AS ic USING(channel_id) '
               'LEFT OUTER JOIN ms1_align USING(spectra_id) '
               'LEFT OUTER JOIN ms1_quant AS pq USING(feature_id) '
               'ORDER BY pr.rownr')
        expected_values = self.process_dbvalues_both(self.dbfile, sql, [1, 2],
                                                     [3], ['MS1 area'])
        fields = ['tmt10plex_{}'.format(ch) for ch in ['126', '127N', '127C',
                                                       '128N', '128C', '129N',
                                                       '129C', '130N', '130C',
                                                       '131']]
        fields.append('MS1 area')
        self.check_results_sql(fields, self.rowify(expected_values))


class TestPercoTSV(basetests.MzidTSVBaseTest):
    command = 'percolator'
    suffix = '_fdr.tsv'

    def test_add_tdc_fdr(self):
        mzidfn = os.path.join(self.fixdir, 'few_spectra.mzid')
        percofn = os.path.join(self.fixdir, 'perco.xml')
        options = ['--mzid', mzidfn, '--perco', percofn]
        self.run_command(options)
        checkfields = ['percolator svm-score', 'PSM q-value', 'peptide q-value', 'TD']
        with open(os.path.join(self.fixdir, 'few_spectra.tsv_fdr.tsv')) as fp:
            header = next(fp).strip().split('\t')
            expected = [line.strip().split('\t') for line in fp]
        expected = [{field: line[i] for i, field in enumerate(header)} for line in expected]
        for res, exp in zip(self.get_values(checkfields), expected):
            for i, field in enumerate(checkfields):
                if field in checkfields:
                    self.assertEqual(field, res[i][1])
                    self.assertEqual(exp[field], res[i][2])


class TestPercoTSVTIMS(basetests.MzidTSVBaseTest):
    command = 'percolator'
    suffix = '_fdr.tsv'
    infilename = 'few_spec_timstof.tsv'
    dbfn = 'spectra_lookup_timstof.sqlite'

    def test_add_tdc_fdr_timstof(self):
        mzidfn = os.path.join(self.fixdir, 'few_spec_timstof.mzid')
        percofn = os.path.join(self.fixdir, 'perco_timstof.xml')
        options = ['--mzid', mzidfn, '--perco', percofn]
        self.run_command(options)
        #checkfields = ['percolator svm-score', 'PSM q-value', 'peptide q-value', 'TD']
        with open(os.path.join(self.fixdir, 'few_spec_timstof.tsv_fdr.tsv')) as fp:
            header = next(fp).strip().split('\t')
            expected = [line.strip().split('\t') for line in fp]
        expected = [{field: line[i] for i, field in enumerate(header)} for line in expected]
        for res, exp in zip(self.get_values(header), expected):
            for i, field in enumerate(header):
                self.assertEqual(field, res[i][1])
                self.assertEqual(exp[field], res[i][2])


class TestMergeTSV(basetests.MzidTSVBaseTest):
    command = 'merge'
    suffix = '_concat.tsv'
    infilename = 'few_spectra.tsv'

    def test_mergetsv(self):
        self.infile = [self.infile, self.infile]
        self.run_command()
        resultlines = self.get_all_lines(self.resultfn)
        for expectfn in self.infile:
            for line in self.get_all_lines(expectfn):
                self.assertEqual(line, next(resultlines))


class TestSplitTSV(basetests.MzidTSVBaseTest):
    infilename = 'few_spectra.tsv'
    command = 'split'
    suffix = '_split.tsv'

    def setUp(self):
        super().setUp()
        self.resultfn = None
        self.expectlines = [x for x in self.get_all_lines(self.infile)]

    def test_auto_bioset_column(self):
        self.run_command(['--bioset'])
        for resultfn in [os.path.join(self.workdir, '{}.tsv'.format(x)) for x in [3, 4, 5]]:
            for line in self.get_all_lines(resultfn):
                self.assertIn(line, self.expectlines)

    def test_splitcol(self):
        setnames = ['3', '4', '5']
        options = ['--splitcol', '8']
        self.run_command(options)
        resultfiles = [os.path.join(self.workdir, '{}.tsv'.format(setname))
                       for setname in setnames]
        for resultfn in resultfiles:
            for line in self.get_all_lines(resultfn):
                self.assertIn(line, self.expectlines)


class TestConffiltTSV(basetests.MzidTSVBaseTest):
    command = 'conffilt'
    infilename = 'few_spectra.tsv'
    suffix = '_filtconf.txt'

    def test_confidence_filter_lower(self):
        conflvl = 0
        self.run_conffilt(conflvl, 'lower', confcol=14)

    def test_confidence_filter_lower_confpattern(self):
        conflvl = 0
        self.run_conffilt(conflvl, 'lower', confpat='EValue')

    def test_confidence_filter_higher(self):
        conflvl = 0
        self.run_conffilt(conflvl, 'higher', confcol=14)

    def run_conffilt(self, conflvl, better, confcol=False, confpat=False):
        options = ['--confidence-better', better,
                   '--confidence-lvl', str(conflvl)]
        if confcol is not False:
            options.extend(['--confidence-col', str(confcol)])
        elif confpat:
            options.extend(['--confcolpattern', confpat])
        self.run_command(options)
        asserter = {'lower': self.assertLess,
                    'higher': self.assertGreater}[better]
        for line in self.get_all_lines(self.resultfn):
            asserter(float(line.strip('\n').split('\t')[confcol - 1]), conflvl)

    def test_confidence_omit_confcol(self):
        options = ['--confidence-better', 'lower', '--confidence-lvl', '0.01']
        self.run_command_expect_error(options)

    def test_omit_conf_better(self):
        options = ['--confidence-col', '1', '--confidence-lvl', '0.01']
        self.run_command_expect_error(options)

    def test_omit_conf_val(self):
        options = ['--confidence-col', '1', '--confidence-better', 'lower']
        self.run_command_expect_error(options)


class TestProteinGroup(basetests.MzidTSVBaseTest):
    command = 'proteingroup'
    infilename = 'target.tsv'
    suffix = '_protgroups.txt'

    def parse_proteingroups(self, fn):
        with open(fn) as fp:
            header = next(fp).strip().split('\t')
            master_ix = header.index(constants.HEADER_MASTER_PROT)
            pgcontent_ix = header.index(constants.HEADER_PG_CONTENT)
            pgamount_ix = header.index(constants.HEADER_PG_AMOUNT_PROTEIN_HITS)
            for line in fp:
                line = line.strip().split('\t')
                yield {'master': line[master_ix],
                       'content': line[pgcontent_ix],
                       'amount': line[pgamount_ix],
                       }

    def do_asserting(self, result, expected, unrolled=False):
        for res, exp in zip(result, expected):
            self.assertEqual(set(res['master'].split(';')),
                             set(exp['master'].split(';')))
            self.assertEqual(res['amount'], exp['amount'])
            rescontent = res['content'].split(';')
            expcontent = exp['content'].split(';')
            self.assertEqual(set(rescontent), set(expcontent))

    def check_lup(self, sql, keyfun, valfun):
        result = {keyfun(x): valfun(x) for x in
                  self.get_values_from_db(self.workdb, sql)}
        exp_file = os.path.join(self.fixdir, 'target_psms_pg.sqlite')
        expected = {keyfun(x): valfun(x) for x in
                    self.get_values_from_db(exp_file, sql)}
        for key, value in result.items():
            self.assertIn(key, expected.keys())
            self.assertEqual(value, expected[key])

    def test_pg(self):
        self.workdb = os.path.join(self.workdir, self.dbfn)
        self.copy_db_to_workdir(self.dbfn, self.workdb)
        options = ['--dbfile', self.workdb]
        self.run_command(options)
        # Check on the lookup
        sql = 'SELECT * FROM protein_coverage'
        self.check_lup(sql, lambda x: x[0], lambda x: x[1])
        sql = ('SELECT ppg.psm_id, pgm.protein_acc FROM psm_protein_groups '
               'AS ppg JOIN protein_group_master AS pgm USING(master_id)')
        self.check_lup(sql, lambda x: x[0], lambda x: x[1])
        sql = ('SELECT pgm.protein_acc, pgc.protein_acc, pgc.peptide_count, '
               'pgc.psm_count, pgc.protein_score '
               'FROM protein_group_content AS pgc '
               'JOIN protein_group_master AS pgm USING(master_id) '
               'ORDER BY pgm.protein_acc, pgc.protein_acc')
        self.check_lup(sql, lambda x: x[0], lambda x: x[1:])
        sql = ('SELECT * FROM protein_group_master')
        self.check_lup(sql, lambda x: x[1], lambda x: 1)
        # Check the output TSV
        result = self.parse_proteingroups(self.resultfn)
        expected = self.parse_proteingroups(
            os.path.join(self.fixdir, 'target_pg.tsv'))
        self.do_asserting(result, expected)


class TestAddGenes(basetests.MzidTSVBaseTest):
    command = 'genes'
    suffix = '_genes.txt'
    infilename = 'target.tsv'
    db_fn = 'target_psms.sqlite'

    def test_addgenes(self):
        self.run_command(['--dbfile', self.dbfile])
        for line in self.get_values(['Gene ID', 'Gene Name', 'Description',
                                     'Protein']):
            genes = line[0][2].split(';')
            assoc_ids = line[1][2].split(';')
            descriptions = ['{}]'.format(x).replace(']]', ']')
                            for x in line[2][2].split('];')]
            proteins = [x.split('(')[0] for x in line[3][2].split(';')]
            sql = ('SELECT p.protein_acc, g.gene_acc, a.assoc_id, '
                   'd.description FROM proteins AS p '
                   'JOIN genes AS g USING(protein_acc) '
                   'JOIN associated_ids AS a ON p.protein_acc=a.protein_acc'
                   ' JOIN prot_desc AS d ON d.protein_acc=p.protein_acc '
                   'WHERE p.protein_acc IN ({})')
            dbvals = self.get_values_from_db(self.dbfile, sql.format(
                ','.join(['"{}"'.format(x) for x in proteins])))
            exp_g, exp_assoc, exp_desc = set(), set(), set()
            for prot, gene, aid, desc in dbvals:
                exp_g.add(gene)
                exp_assoc.add(aid)
                exp_desc.add(desc)
            for exp_set, result in zip([exp_g, exp_assoc, exp_desc],
                                       [genes, assoc_ids, descriptions]):
                self.assertEqual(0, len(exp_set.difference(result)))


class TestIso(basetests.MzidTSVBaseTest):

    def get_denominator(self, line, denom_ch):
        denomvals = [float(line[ch]) for ch in denom_ch if line[ch] != 'NA']
        return sum(denomvals) / len(denomvals)

    def get_infile_lines(self, infile=None):
        if infile is None:
            infile = self.infile[0]
        with open(infile) as fp:
            header = next(fp).strip('\n').split('\t')
            for line in fp:
                line = line.strip('\n').split('\t')
                yield {field: val for field, val in zip(header, line)}

    def check_normalize_medians(self, channels, denom_ch, minint, stdout,
                                medianpsms):
        ch_medians = {ch: [] for ch in channels}
        for line in self.get_infile_lines(medianpsms):
            line.update({ch: line[ch]
                         if line[ch] != 'NA' and float(line[ch]) > minint
                         else 'NA' for ch in channels})
            denom = self.get_denominator(line, denom_ch)
            if denom == 0:
                continue
            for ch in channels:
                if line[ch] == 'NA':
                    continue
                ch_medians[ch].append(float(line[ch]) / denom)
        ch_medians = {ch: median(vals) for ch, vals in ch_medians.items()}
        stdout = stdout.decode().split('\n')
        self.assertEqual(stdout[0],
                         'Channel intensity medians used for normalization:')
        stdout_channels = {x.split(' - ')[0]: x.split(' - ')[1]
                           for x in stdout[1:]}
        for ch in channels:
            self.assertEqual(float(stdout_channels[ch]), ch_medians[ch])
        return ch_medians

    def do_check(self, minint, stdout, normalize=False, medianpsms=None,
                 resultch=False):
        channels = ['tmt10plex_126'] + [x.format('tmt10plex_1', y+27) for x in ['{}{}C', '{}{}N'] for y in range(4)] + ['tmt10plex_131']
        resultch = ['ratio_{}'.format(x) for x in channels]
        denom_ch = [channels[0], channels[-1]]
        if normalize:
            ch_medians = self.check_normalize_medians(channels, denom_ch,
                                                      minint, stdout,
                                                      medianpsms)
        for in_line, resultline in zip(self.get_infile_lines(),
                                       self.get_values(resultch)):
            in_line.update({ch: in_line[ch]
                            if in_line[ch] != 'NA' and
                            float(in_line[ch]) > minint else 'NA'
                            for ch in channels})
            resultline = [x[2] for x in resultline]
            denom = self.get_denominator(in_line, denom_ch)
            if denom == 0:
                exp_line = ['NA'] * len(channels)
            elif normalize:
                exp_line = [str((float(in_line[ch]) / denom) / ch_medians[ch])
                            if in_line[ch] != 'NA' else 'NA'
                            for ch in channels]
            else:
                exp_line = [str((float(in_line[ch]) / denom))
                            if in_line[ch] != 'NA' else 'NA'
                            for ch in channels]
            self.assertEqual(resultline, exp_line)


class TestIsoRatio(TestIso):
    suffix = '_ratio_isobaric.txt'
    command = 'isoratio'
    infilename = 'quant_target.tsv'

    def test_denomcolpattern(self):
        stdout = self.run_command_stdout(['--isobquantcolpattern', 'plex',
                                          '--denompatterns', '_126', '_131'])
        self.do_check(0, stdout)

    def test_denomcolpattern_regex(self):
        stdout = self.run_command_stdout(['--isobquantcolpattern', 'plex',
                                          '--denompatterns', '_1[23][61]'])
        self.do_check(0, stdout)


class TestIsoFeatRatio(TestIso):
    suffix = '_ratio_isobaric.txt'
    command = 'isoratio'
    infilename = 'target_pg.tsv'
    channels = ['tmt10plex_{}'.format(x) for x in ['126', '127N', '127C',
                                                   '128N', '128C', '129N',
                                                   '129C', '130N', '130C',
                                                   '131']]
    nopsms = ['{} - # quanted PSMs'.format(ch) for ch in channels]

    def test_normalized_isoquant(self):
        options = ['--protcol', '11', '--isobquantcolpattern', 'tmt10plex',
                   '--denompatterns', '_126', '--normalize', 'median']
        self.run_command(options)
        self.isoquant_check(os.path.join(self.fixdir, 'proteins_normquant'),
            'Protein ID', self.channels, self.nopsms)

    def test_normalized_targettable_isoquant(self):
        prottable_targettable = os.path.join(self.fixdir, 'target_proteins')
        options = ['--protcol', '11', '--isobquantcolpattern', 'tmt10plex',
                   '--denompatterns', '_126', '--normalize', 'median',
                   '--targettable', prottable_targettable]
        self.run_command(options)
        self.isoquant_check(os.path.join(self.fixdir, 'proteins_normquant'),
                'Protein ID', self.channels, self.nopsms)
