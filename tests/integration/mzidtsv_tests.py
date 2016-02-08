import os
from lxml import etree

from app.dataformats import mzidtsv as constants
from tests.integration import basetests


class TestAddSpecData(basetests.MzidTSVBaseTest):
    command = 'addspecdata'
    suffix = '_spectradata.tsv'
    infilename = 'mzidtsv_filtered_fr1-2_nospecdata.txt'

    def test_addspecdata(self):
        options = ['--dbfile', self.dbfile, '--spectracol', '2']
        self.run_command(options)
        sql = ('SELECT pr.rownr, bs.set_name, sp.retention_time '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'JOIN mzml AS sp USING(spectra_id) '
               'JOIN mzmlfiles USING(mzmlfile_id) '
               'JOIN biosets AS bs USING(set_id) '
               'ORDER BY pr.rownr')
        fields = ['Biological set', 'Retention time(min)']
        expected_values = self.process_dbvalues_both(self.dbfile, sql, [],
                                                     [1, 2], fields)
        self.check_results(fields, expected_values)


class TestQuantTSV(basetests.MzidTSVBaseTest):
    command = 'addquant'
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
        self.check_results(fields, expected_values)

    def test_quanttsv_precursor(self):
        dbfile = os.path.join(self.fixdir, 'mzidtsv_db.sqlite')
        options = ['--dbfile', dbfile, '--precursor']
        self.run_command(options)
        sql = ('SELECT pr.rownr, pq.intensity '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'LEFT OUTER JOIN ms1_align USING(spectra_id) '
               'LEFT OUTER JOIN ms1_quant AS pq USING(feature_id)')
        expected_values = self.get_values_from_db(self.dbfile, sql)
        self.check_results(['MS1 area'], expected_values)

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
        expected_values = self.process_dbvalues_both(self.dbfile, sql,
                                                     [1, 2], [3], ['MS1 area'])
        fields = ['tmt10plex_{}'.format(ch) for ch in ['126', '127N', '127C',
                                                       '128N', '128C', '129N',
                                                       '129C', '130N', '130C',
                                                       '131']]
        fields.append('MS1 area')
        self.check_results(fields, expected_values)


class TestPercoTSV(basetests.MzidTSVBaseTest):
    command = 'addperco'
    suffix = '_percolated.tsv'
    infilename = 'mzidtsv_fr0.txt'
    field_p_map = {'percolator svm-score': 'score',
                   'PSM p-value': 'psm_p_value',
                   'PSM q-value': 'psm_q_value',
                   'PSM-PEP': 'psm_pep',
                   'peptide q-value': 'peptide_q_value',
                   'peptide PEP': 'peptide PEP',
                   }

    def test_add_percolator(self):
        mzidfn = os.path.join(self.fixdir, 'msgfperco_fr0.mzid')
        options = ['--mzid', mzidfn]
        self.run_command(options)
        expected = self.get_percolator_from_msgf(mzidfn,
                                                 self.field_p_map.keys())
        self.check_results(self.field_p_map.keys(), expected)

    def get_percolator_from_msgf(self, msgffile, checkfields):
        for specidres in etree.iterparse(msgffile,
                                         tag='SpectrumIdentificationResult'):
            count = 0
            for result in specidres:
                sid = result.attrib['spectrumID']
                perco = [x for x in result.findall('userParam')
                         if x.attrib['name'].split(':')[0] == 'percolator']
                perco = {x.attrib['name'].replace('percolator:', ''):
                         x.attrib['value'] for x in perco}
                outresult = [(count, 'SpecID', sid)]
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
        self.expectlines = self.get_all_lines(self.infile)

    def test_auto_bioset_column(self):
        self.run_command(['--bioset', '--setnames', 'S1'])
        print(os.listdir(self.workdir))
        resultfn = os.path.join(self.workdir, '{}_0{}'.format(
            os.path.basename(self.infile[0]), self.suffix))
        for line in self.get_all_lines(resultfn):
            self.assertEqual(line, next(self.expectlines))

    def test_splitcol(self):
        self.run_command(['--splitcol', '1', '--setnames',
                          'dataset_17694.dat_task_0.mzml',
                          'dataset_17694.dat_task_1.mzml'])
        resultfiles = [os.path.join(self.workdir,
                                    '{}_{}{}'.format(self.infilename,
                                                     index, self.suffix))
                       for index in [0, 1]]
        for resultfn in resultfiles:
            for line in self.get_all_lines(resultfn):
                self.assertEqual(line, next(self.expectlines))


class TestConffiltTSV(basetests.MzidTSVBaseTest):
    command = 'conffilt'
    infilename = 'mzidtsv_fr0.txt'
    suffix = '_filtconf.txt'

    def test_confidence_filter_lower(self):
        conflvl, confcol = 0, 14  # EValue
        self.run_conffilt(confcol, conflvl, True)

    def test_confidence_filter_higher(self):
        conflvl, confcol = 0, 14  # EValue
        self.run_conffilt(confcol, conflvl, False)

    def run_conffilt(self, confcol, conflvl, lowerbetter):
        options = ['--confidence-col', str(confcol), '--confidence-better',
                   'lower', '--confidence-lvl', str(conflvl)]
        self.run_command(options)
        asserter = {True: self.assertLess,
                    False: self.assertGreater}[lowerbetter]
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
