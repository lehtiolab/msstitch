import os
from lxml import etree
from statistics import median

from app.dataformats import mzidtsv as constants
from tests.integration import basetests


class TestAddSpecData(basetests.MzidTSVBaseTest):
    command = 'specdata'
    suffix = '_spectradata.tsv'
    infilename = 'mzidtsv_filtered_fr1-2_nospecdata.txt'

    def test_addspecdata(self):
        options = ['--dbfile', self.dbfile, '--spectracol', '2']
        self.run_command(options)
        sql = ('SELECT pr.rownr, bs.set_name, sp.retention_time, '
               'sp.ion_injection_time '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'JOIN mzml AS sp USING(spectra_id) '
               'JOIN mzmlfiles USING(mzmlfile_id) '
               'JOIN biosets AS bs USING(set_id) '
               'ORDER BY pr.rownr')
        fields = ['Biological set', 'Retention time(min)',
                  'Ion injection time(ms)']
        expected_values = self.process_dbvalues_both(self.dbfile, sql, [],
                                                     [1, 2, 3], fields)
        self.check_results(fields, self.rowify(expected_values))


class TestQuantTSV(basetests.MzidTSVBaseTest):
    command = 'quant'
    suffix = '_quant.tsv'

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
        self.check_results(fields, self.rowify(expected_values))

    def test_quanttsv_precursor(self):
        dbfile = os.path.join(self.fixdir, 'mzidtsv_db.sqlite')
        options = ['--dbfile', dbfile, '--precursor']
        self.run_command(options)
        sql = ('SELECT pr.rownr, pq.intensity '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'LEFT OUTER JOIN ms1_align USING(spectra_id) '
               'LEFT OUTER JOIN ms1_quant AS pq USING(feature_id)')
        expected_values = self.get_values_from_db(self.dbfile, sql)
        self.check_results(['MS1 area'], self.rowify(expected_values))

    def test_quanttsv_both(self):
        dbfile = os.path.join(self.fixdir, 'mzidtsv_db.sqlite')
        options = ['--dbfile', dbfile, '--isobaric', '--precursor']
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
        self.check_results(fields, self.rowify(expected_values))


class TestPercoTSV(basetests.MzidTSVBaseTest):
    command = 'percolator'
    suffix = '_percolated.tsv'
    infilename = 'mzidtsv_fr0.txt'
    field_p_map = {'percolator svm-score': 'score',
                   'PSM p-value': 'psm_p_value',
                   'PSM q-value': 'psm_q_value',
                   'PSM-PEP': 'psm_pep',
                   'peptide q-value': 'peptide_q_value',
                   'peptide PEP': 'peptide_pep',
                   }

    def test_add_percolator(self):
        mzidfn = os.path.join(self.fixdir, 'msgfperco_fr0.mzid')
        options = ['--mzid', mzidfn]
        self.run_command(options)
        expected = self.get_percolator_from_msgf(mzidfn,
                                                 self.field_p_map.keys())
        self.check_results(self.field_p_map.keys(), expected)

    def get_percolator_from_msgf(self, msgffile, checkfields):
        count = 0
        ns = self.get_xml_namespace(msgffile)
        for ac, specidres in etree.iterparse(
                msgffile, tag='{%s}'
                'SpectrumIdentificationResult' % ns['xmlns']):
            for result in specidres.findall(
                    '{%s}SpectrumIdentificationItem' % ns['xmlns']):
                perco = [x for x in
                         result.findall('{%s}userParam' % ns['xmlns'])
                         if x.attrib['name'].split(':')[0] == 'percolator']
                perco = {x.attrib['name'].replace('percolator:', ''):
                         x.attrib['value'] for x in perco}
                if not perco:
                    perco = {key: None for key in self.field_p_map.values()}
                #outresult = [(count, 'SpecID', sid)]
                outresult = []
                outresult.extend([(count, field,
                                   perco[self.field_p_map[field]])
                                  for field in checkfields])
                yield outresult
                count += 1


class TestMergeTSV(basetests.MzidTSVBaseTest):
    command = 'merge'
    suffix = '_concat.tsv'
    infilename = 'mzidtsv_fr0.txt'

    def test_mergetsv(self):
        self.infile = [self.infile, os.path.join(self.fixdir,
                                                 'mzidtsv_fr1.txt')]
        self.run_command()
        resultlines = self.get_all_lines(self.resultfn)
        for expectfn in self.infile:
            for line in self.get_all_lines(expectfn):
                self.assertEqual(line, next(resultlines))


class TestSplitTSV(basetests.MzidTSVBaseTest):
    infilename = 'mzidtsv_filtered_fr1-2.txt'
    command = 'split'
    suffix = '_split.tsv'

    def setUp(self):
        super().setUp()
        self.resultfn = None
        self.expectlines = self.get_all_lines(self.infile)

    def test_auto_bioset_column(self):
        self.run_command(['--bioset'])
        resultfn = os.path.join(self.workdir, 'S1.tsv')
        for line in self.get_all_lines(resultfn):
            self.assertEqual(line, next(self.expectlines))

    def test_splitcol(self):
        setnames = ['dataset_17694.dat_task_0.mzml',
                    'dataset_17694.dat_task_1.mzml']
        options = ['--splitcol', '1']
        self.run_command(options)
        resultfiles = [os.path.join(self.workdir, '{}.tsv'.format(setname))
                       for setname in setnames]
        for resultfn in resultfiles:
            for line in self.get_all_lines(resultfn):
                self.assertEqual(line, next(self.expectlines))


class TestConffiltTSV(basetests.MzidTSVBaseTest):
    command = 'conffilt'
    infilename = 'mzidtsv_fr0.txt'
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
        print(options)
        self.run_command(options)
        asserter = {'lower': self.assertLess,
                    'higher': self.assertGreater}[better]
        for line in self.get_all_lines(self.resultfn):
            asserter(float(line.strip('\n').split('\t')[confcol - 1]), conflvl)

    def test_confidence_omit_confcol(self):
        options = ['--confidence-better', 'lower', '--confidence-lvl', '0.01']
        self.run_command_expect_error(options)

    def test_omit_conf_better(self):
        options = ['--confidence-col', '18', '--confidence-lvl', '0.01']
        self.run_command_expect_error(options)

    def test_omit_conf_val(self):
        options = ['--confidence-col', '18', '--confidence-better', 'lower']
        self.run_command_expect_error(options)


class TestProteinGroup(basetests.MzidTSVBaseTest):
    command = 'proteingroup'
    infilename = 'mzidtsv_filtered_fr1-2.txt'
    suffix = '_protgroups.txt'
    dbfile = 'mzidtsv_db.sqlite'

    def run_and_analyze(self, options):
        self.run_command(options)
        result = self.parse_proteingroups(self.resultfn)
        expected = self.parse_proteingroups(
            os.path.join(self.fixdir,
                         'mzidtsv_filtered_fr1-2_proteingrouped.txt'))
        self.do_asserting(result, expected)

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

    def test_proteingroups(self):
        options = ['--dbfile', self.dbfile]
        self.expected = None
        self.run_and_analyze(options)


class TestAddGenes(basetests.MzidTSVBaseTest):
    command = 'genes'
    suffix = '_genes.txt'
    infilename = 'mzidtsv_filtered_fr1-2.txt'

    def test_addgenes(self):
        self.run_command(['--dbfile', self.dbfile])
        for line in self.get_values(['Gene', 'Gene Symbol', 'Description',
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
        channels = ['fake_ch{}'.format(x) for x in range(8)]
        # TODO only for backwards compatibilty, remove if statement around
        # assignment when msspsmtable isonormalize is removed
        if not resultch:
            resultch = ['ratio_{}'.format(x) for x in channels]
        denom_ch = channels[0:2]
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
    infilename = 'mzidtsv.txt'

    def test_denomcolpattern(self):
        stdout = self.run_command_stdout(['--isobquantcolpattern', 'fake_ch',
                                          '--denompatterns', '_ch0', '_ch1'])
        self.do_check(0, stdout)

    def test_denomcolpattern_regex(self):
        stdout = self.run_command_stdout(['--isobquantcolpattern', 'fake_ch',
                                          '--denompatterns', '_ch[0-1]'])
        self.do_check(0, stdout)


class TestIsoFeatRatio(TestIso):
    suffix = '_ratio_isobaric.txt'
    command = 'isoratio'
    infilename = 'mzidtsv_intensities.txt'
    channels = ['tmt10plex_{}'.format(x) for x in ['126', '127N', '127C',
                                                   '128N', '128C', '129N',
                                                   '129C', '130N', '130C',
                                                   '131']]
    nopsms = ['{} - # quanted PSMs'.format(ch) for ch in channels]

    def test_normalized_isoquant(self):
        options = ['--protcol', '14', '--isobquantcolpattern', 'tmt10plex',
                   '--denompatterns', '_126', '--normalize', 'median']
        self.run_command(options)
        self.isoquant_check(
            os.path.join(self.fixdir, 'prottable_normalized_isoquant.txt'),
            'Accession', self.channels, self.nopsms)

    def test_normalized_othertable_isoquant(self):
        prottable_ratiofn = os.path.join(self.fixdir,
                                         'mzidtsv_ratios.txt')
        options = ['--protcol', '14', '--isobquantcolpattern', 'tmt10plex',
                   '--denompatterns', '_126', '--normalize', 'median',
                   '--norm-ratios', prottable_ratiofn]
        self.run_command(options)
        self.isoquant_check(
            os.path.join(self.fixdir, 'prottable_normalized_isoquant.txt'),
            'Accession', self.channels, self.nopsms)

    def test_normalized_targettable_isoquant(self):
        prottable_targettable = os.path.join(self.fixdir,
                                             'prottable_only_acc.txt')
        options = ['--protcol', '14', '--isobquantcolpattern', 'tmt10plex',
                   '--denompatterns', '_126', '--normalize', 'median',
                   '--targettable', prottable_targettable]
        self.run_command(options)
        self.isoquant_check(
            os.path.join(self.fixdir, 'prottable_normalized_isoquant.txt'),
            'Protein accession', self.channels, self.nopsms)


class TestIsoNormalize(TestIso):
    suffix = '_normalized_isobaric.txt'
    command = 'isonormalize'
    infilename = 'mzidtsv.txt'
    channels = ['fake_ch{}'.format(x) for x in range(8)]

    def test_normalize(self):
        stdout = self.run_command_stdout(['--isobquantcolpattern', 'fake_ch',
                                          '--denomcols', '21', '22'])
        self.do_check(0, stdout, normalize=True, resultch=self.channels)

    def test_normalize_minint(self):
        minint = 3000
        stdout = self.run_command_stdout(['--isobquantcolpattern', 'fake_ch',
                                          '--denomcols', '21', '22',
                                          '--minint', str(minint)])
        self.do_check(minint, stdout, normalize=True, resultch=self.channels)


class TestIsoNormalizeTwofiles(TestIso):
    infilename = 'mzidtsv_short.txt'
    suffix = '_normalized_isobaric.txt'
    command = 'isonormalize'
    channels = ['fake_ch{}'.format(x) for x in range(8)]

    def test_two_psm_files(self):
        """Tests calculating medians on different file than the one doing the
        median centering on"""
        medianpsms = os.path.join(self.fixdir, 'mzidtsv.txt')
        stdout = self.run_command_stdout(['--isobquantcolpattern', 'fake_ch',
                                          '--denomcols', '21', '22',
                                          '--medianpsms', medianpsms])
        self.do_check(0, stdout, normalize=True, medianpsms=medianpsms,
                      resultch=self.channels)
