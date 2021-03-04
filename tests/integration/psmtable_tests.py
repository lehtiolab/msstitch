import os
import re
import subprocess
from lxml import etree
from Bio import SeqIO
from statistics import median

from app.dataformats import mzidtsv as constants
from tests.integration import basetests


class MzidWithDB(basetests.MzidTSVBaseTest):
    def setUp(self):
        super().setUp()
        self.dbfile = os.path.join(self.fixdir, self.dbfn)
        self.workdb = os.path.join(self.workdir, self.dbfn)
        self.copy_db_to_workdir(self.dbfn, self.workdb)


class TestPSMTable(MzidWithDB):
    command = 'psmtable'
    infilename = 'target.tsv'
    dbfn = 'quant_lookup.sqlite'
    """DB and PSM table contain:
    - ENSEMBL proteins
    - a Uniprot swiss protein
    - A self-annotated protein
    - A non-annotated (only peptide) proteins
    """

    def test_build_full_psmtable(self):
        minpif = '0.4'
        fastafn = os.path.join(self.basefixdir, 'ens99_small.fasta')
        options = ['--dbfile', self.workdb, '--spectracol', '1', '--addmiscleav',
                '--addbioset', '--genes', '--proteingroup', '--ms1quant', '--isobaric',
                '--fasta', fastafn, '--min-precursor-purity', minpif]
        self.run_command(options)
        self.check_db_fasta(fastafn)
        self.check_addspec_miscleav_bioset()
        self.check_pg()
        self.check_quanttsv(minpif)
        self.check_addgenes()

    def test_ionmobility(self):
        self.infilename = 'few_spec_timstof.tsv'
        self.dbfn = 'spectra_lookup_timstof.sqlite'
        self.setUp()
        options = ['--dbfile', self.workdb, '--spectracol', '1', '--addmiscleav', '--addbioset']
        self.run_command(options)
        sql = ('SELECT pr.rownr, bs.set_name, sp.retention_time, '
               'iit.ion_injection_time, im.ion_mobility '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'JOIN mzml AS sp USING(spectra_id) '
               'LEFT OUTER JOIN ioninjtime AS iit USING(spectra_id) '
               'LEFT OUTER JOIN ionmob AS im USING(spectra_id) '
               'JOIN mzmlfiles USING(mzmlfile_id) '
               'JOIN biosets AS bs USING(set_id) '
               'ORDER BY pr.rownr')
        fields = ['Biological set', 'Retention time(min)',
                  'Ion injection time(ms)', 'Ion mobility(Vs/cm2)']
        expected_values = self.process_dbvalues_both(self.workdb, sql, [],
                                                     [1, 2, 3, 4], fields)
        self.check_results_sql(fields, self.rowify(expected_values))
        for val, exp in zip(self.get_values(['missed_cleavage']), self.get_values(['Peptide'], self.infile[0])):
            exp = re.sub('[0-9\+\.]', '', exp[0][1])[:-1]
            self.assertEqual(int(val[0][1]), exp.count('K') + exp.count('R') - exp.count('KP') - exp.count('RP'))

    def check_addspec_miscleav_bioset(self):
        sql = ('SELECT pr.rownr, bs.set_name, sp.retention_time, '
               'iit.ion_injection_time, im.ion_mobility, pif.pif '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'JOIN mzml AS sp USING(spectra_id) '
               'JOIN ioninjtime AS iit USING(spectra_id) '
               'LEFT OUTER JOIN ionmob AS im USING(spectra_id) '
               'LEFT OUTER JOIN precursor_ion_fraction AS pif USING(spectra_id) '
               'JOIN mzmlfiles USING(mzmlfile_id) '
               'JOIN biosets AS bs USING(set_id) '
               'ORDER BY pr.rownr')
        fields = ['Biological set', 'Retention time(min)',
                  'Ion injection time(ms)', 'Ion mobility(Vs/cm2)',
                  'Precursor ion fraction']
        expected_values = self.process_dbvalues_both(self.workdb, sql, [],
                                                     [1, 2, 3, 4, 5], fields)
        self.check_results_sql(fields, self.rowify(expected_values))
        for val, exp in zip(self.get_values(['missed_cleavage']), self.get_values(['Peptide'], self.infile[0])):
            exp = re.sub('[0-9\+\.]', '', exp[0][1])[:-1]
            self.assertEqual(int(val[0][1]), exp.count('K') + exp.count('R') - exp.count('KP') - exp.count('RP'))

    def check_quanttsv(self, minpif):
        sql = ('SELECT pr.rownr, ic.channel_name, '
                'CASE WHEN pif.pif > {} THEN iq.intensity ELSE "NA" END AS intensity '
                'FROM psmrows AS pr JOIN psms USING(psm_id) '
                'LEFT OUTER JOIN precursor_ion_fraction AS pif USING(spectra_id) '
                'JOIN isobaric_quant AS iq USING(spectra_id) '
                'JOIN isobaric_channels AS ic USING(channel_id) '.format(float(minpif)))
        expected_values = self.get_values_from_db(self.workdb, sql)
        fields = ['tmt10plex_{}'.format(ch) for ch in ['126', '127N', '127C',
                                                       '128N', '128C', '129N',
                                                       '129C', '130N', '130C',
                                                       '131']]
        self.check_results_sql(fields, self.rowify(expected_values))
        sql = ('SELECT pr.rownr, pq.intensity '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'LEFT OUTER JOIN ms1_align USING(spectra_id) '
               'LEFT OUTER JOIN ms1_quant AS pq USING(feature_id)')
        expected_values = self.get_values_from_db(self.workdb, sql)
        self.check_results_sql(['MS1 area'], self.rowify(expected_values))

    def check_db_fasta(self, fasta, exp_proteins=None, desc=True):
        if exp_proteins is None:
            exp_proteins = {}
            for rec in SeqIO.parse(fasta, 'fasta'):
                rd = rec.description
                gene = 'NA'
                if 'gene_symbol:' in rd:
                    six = rd.index('gene_symbol:') + 12
                    gix = rd.index('gene:') + 5 
                    gene = rd[gix: rd.index(' ', gix)]
                    desc = rd[rd.index('description:') + 12:]
                elif 'GN=' in rd:
                    six = rd.index('GN=') + 3
                    desc = [x for x in rd.split() if '=' not in x][1:]
                elif 'msstitch_fake_gene' in rd:
                    # special case fake fasta record for non-standard gene 
                    six, desc = False, 'NA'
                    symbol = rd.split()[-1]
                elif 'msstitch_fake_onlypeptide' in rd:
                    # special fake fasta record for unannotated peptide
                    six, symbol, desc = False, 'NA', 'NA'
                exp_proteins[rec.id] = {
                        'seq': rec.seq,
                        'gene': gene,
                        'desc': desc,
                        'symbol': rd[six: rd.index(' ', six)] if six else 'NA',
                        }
        self.check_db_base(exp_proteins)
        sql = ('SELECT ps.protein_acc, ps.sequence, g.gene_acc, aid.assoc_id, '
                'pd.description '
                'FROM genes AS g '
                'JOIN ensg_proteins USING(gene_id) '
                'JOIN genename_proteins AS gnp USING(pacc_id) '
                'JOIN associated_ids AS aid USING(gn_id) '
                'JOIN proteins AS p USING(pacc_id) '
               ' JOIN protein_seq AS ps ON ps.protein_acc=p.protein_acc '
               ' JOIN prot_desc AS pd ON pd.pacc_id=p.pacc_id ')
        if not desc:
            sql = ('SELECT ps.protein_acc, ps.sequence '
                   'FROM protein_seq AS ps '
                   'JOIN prot_desc AS pd USING(protein_acc)')
        for prot, seq, gene, aid, desc in self.get_values_from_db(self.workdb,
                                                             sql):
            self.assertEqual(exp_proteins[prot]['seq'], seq)
            if desc:
                self.assertEqual(exp_proteins[prot]['gene'], gene)
                self.assertEqual(exp_proteins[prot]['symbol'], aid)
                self.assertEqual(exp_proteins[prot]['desc'], desc)

    def check_db_base(self, expected_proteins=None):
        expected_psms = self.get_expected_psms()
        if expected_proteins is None:
            expected_proteins = {x for y in expected_psms.values()
                                 for x in y['proteins']}
        protsql = 'SELECT * FROM proteins'
        for protein in self.get_values_from_db(self.workdb, protsql):
            self.assertIn(protein[1], expected_proteins)
        psmsql = ('SELECT ps.sequence, p.score, pr.rownr '
                  'FROM psmrows AS pr JOIN psms AS p USING(psm_id) '
                  'JOIN peptide_sequences AS ps USING(pep_id)')
        for psm in self.get_values_from_db(self.workdb, psmsql):
            expected_psm = (expected_psms[psm[2]]['seq'],
                            expected_psms[psm[2]]['score'])
            self.assertEqual((psm[0], psm[1]), expected_psm)
        ppsql = ('SELECT pp.protein_acc, pr.rownr FROM psmrows AS pr '
                 'JOIN protein_psm AS pp USING(psm_id)')
        for protpsm in self.get_values_from_db(self.workdb, ppsql):
            self.assertIn(protpsm[0], expected_psms[protpsm[1]]['proteins'])

    def get_expected_psms(self):
        header = self.get_tsvheader(self.infile[0])
        prot_ix = header.index('Protein')
        seq_ix = header.index('Peptide')
        score_ix = header.index('MSGFScore')
        psms = {}
        for row, line in enumerate(self.get_all_lines(self.infile[0])):
            line = line.strip('\n').split('\t')
            psms[row] = {'proteins': [x.split('(pre')[0] for x in
                                      line[prot_ix].split(';')],
                         'seq': line[seq_ix],
                         'score': line[score_ix],
                         }
        return psms

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

    def check_pglup(self, sql, keyfun, valfun):
        result = {keyfun(x): valfun(x) for x in
                  self.get_values_from_db(self.workdb, sql)}
        exp_file = os.path.join(self.fixdir, 'target_psms.sqlite')
        expected = {keyfun(x): valfun(x) for x in
                    self.get_values_from_db(exp_file, sql)}
        for key, value in result.items():
            self.assertIn(key, expected.keys())
            self.assertEqual(value, expected[key])

    def check_pg(self):
        sql = 'SELECT * FROM protein_coverage'
        self.check_pglup(sql, lambda x: x[0], lambda x: x[1])
        sql = """SELECT ppg.psm_id, p.protein_acc FROM psm_protein_groups
               AS ppg JOIN protein_group_master AS pgm USING(master_id)
               JOIN proteins AS p ON pgm.pacc_id=p.pacc_id"""
        self.check_pglup(sql, lambda x: x[0], lambda x: x[1])
        sql = ('SELECT p.protein_acc, pgc.protein_acc, pgc.peptide_count, '
               'pgc.psm_count, pgc.protein_score '
               'FROM protein_group_content AS pgc '
               'JOIN protein_group_master AS pgm USING(master_id) '
               'JOIN proteins AS p USING(pacc_id) '
               'ORDER BY p.pacc_id, pgc.protein_acc')
        self.check_pglup(sql, lambda x: x[0], lambda x: x[1:])
        sql = ('SELECT * FROM protein_group_master')
        self.check_pglup(sql, lambda x: x[1], lambda x: 1)
        # Check the output TSV
        result = self.parse_proteingroups(self.resultfn)
        expected = self.parse_proteingroups(
            os.path.join(self.fixdir, 'target_pg.tsv'))
        self.do_asserting(result, expected)

    def check_addgenes(self):
        for line in self.get_values(['Gene ID', 'Gene Name', 'Description',
                                     'Protein']):
            genes = line[0][2].split(';')
            assoc_ids = line[1][2].split(';')
            descriptions = ['{}]'.format(x).replace(']]', ']')
                            for x in line[2][2].split('];')]
            proteins = [x.split('(')[0] for x in line[3][2].split(';')]
            sql = ('SELECT p.protein_acc, g.gene_acc, a.assoc_id, '
                   'd.description FROM proteins AS p '
                   'JOIN ensg_proteins USING(pacc_id) '
                   'JOIN genename_proteins USING(pacc_id) '
                   'JOIN genes AS g USING(gene_id) '
                   'JOIN associated_ids AS a USING(gn_id) '
                   ' JOIN prot_desc AS d ON d.pacc_id=p.pacc_id '
                   'WHERE p.protein_acc IN ({})')
            dbvals = self.get_values_from_db(self.workdb, sql.format(
                ','.join(['"{}"'.format(x) for x in proteins])))
            exp_g, exp_assoc, exp_desc = set(), set(), set()
            for prot, gene, aid, desc in dbvals:
                exp_g.add(gene)
                exp_assoc.add(aid)
                exp_desc.add(desc)
            for exp_set, result in zip([exp_g, exp_assoc, exp_desc],
                                       [genes, assoc_ids, descriptions]):
                self.assertEqual(0, len(exp_set.difference(result)))


class TestPercoTSV(basetests.MzidTSVBaseTest):
    command = 'perco2psm'
    suffix = '_fdr.tsv'
    infilename = 'few_spectra.tsv'

    def test_conffilt(self):
        mzidfn = os.path.join(self.fixdir, 'few_spectra.mzid')
        percofn = os.path.join(self.fixdir, 'perco.xml')
        options = ['--mzid', mzidfn, '--perco', percofn, '--filtpep', '0.01',
            '--filtpsm', '0.01']
        self.run_command(options)
        checkfields = ['percolator svm-score', 'PSM q-value', 'peptide q-value', 'TD']
        with open(os.path.join(self.fixdir, 'few_spectra.tsv_fdr.tsv')) as fp:
            header = next(fp).strip().split('\t')
            expected = [line.strip().split('\t') for line in fp]
        expected = [{field: line[i] for i, field in enumerate(header)} for line in expected]
        expected = [line for line in expected if float(line['PSM q-value']) < 0.01 and
            float(line['peptide q-value']) < 0.01]
        self.assertEqual(len(expected),  len([x for x in self.get_values(checkfields)]))
        for res, exp in zip(self.get_values(checkfields), expected):
            for i, field in enumerate(checkfields):
                if field in checkfields:
                    self.assertEqual(field, res[i][1])
                    self.assertEqual(exp[field], res[i][2])

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
    command = 'perco2psm'
    suffix = '_fdr.tsv'
    infilename = 'few_spec_timstof.tsv'
    dbfn = 'spectra_lookup_timstof.sqlite'

    def test_add_tdc_fdr_timstof(self):
        mzidfn = os.path.join(self.fixdir, 'few_spec_timstof.mzid')
        percofn = os.path.join(self.fixdir, 'perco_timstof.xml')
        options = ['--mzid', mzidfn, '--perco', percofn]
        self.run_command(options)
        with open(os.path.join(self.fixdir, 'few_spec_timstof.tsv_fdr.tsv')) as fp:
            header = next(fp).strip().split('\t')
            expected = [line.strip().split('\t') for line in fp]
        expected = [{field: line[i] for i, field in enumerate(header)} for line in expected]
        for res, exp in zip(self.get_values(header), expected):
            for i, field in enumerate(header):
                self.assertEqual(field, res[i][1])
                self.assertEqual(exp[field], res[i][2])


class TestConcatTSV(basetests.MzidTSVBaseTest):
    command = 'concat'
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
    infilename = 'target_pg.tsv'
    command = 'split'
    suffix = '_split.tsv'

    def setUp(self):
        super().setUp()
        self.resultfn = None
        self.expectlines = [x for x in self.get_all_lines(self.infile)]

    def test_splitcol_bioset(self):
        self.run_command(['--splitcol', 'bioset'])
        for resultfn in [os.path.join(self.workdir, '{}.tsv'.format(x)) for x in ['Set1', 'Set2']]:
            for line in self.get_all_lines(resultfn):
                self.assertIn(line, self.expectlines)

    def test_invalid_splitcol(self):
        options = ['--splitcol', 'hello']
        result = self.run_command(options, return_error=True)
        if result.returncode != 0:
            self.assertEqual(result.stdout.strip(), 
                    'ERROR: --splitcol must be an integer or "TD", or "bioset"')
        else:
            self.fail('This test should error')

    def test_splitcol(self):
        setnames = ['Set1', 'Set2']
        options = ['--splitcol', '28']
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
        res = self.run_command(options, return_error=True)
        if res.returncode != 0:
            self.assertEqual(res.stdout.strip(), 
                    'Must define either --confcol or --confcolpattern')
        else:
            self.fail('This test should error')

    def test_omit_conf_better(self):
        options = ['--confidence-col', '1', '--confidence-lvl', '0.01']
        res = self.run_command(options, return_error=True)
        if res.returncode != 0:
            self.assertIn('usage', res.stderr.strip())
        else:
            self.fail('This test should error')

    def test_omit_conf_val(self):
        options = ['--confidence-col', '1', '--confidence-better', 'lower']
        res = self.run_command(options, return_error=True)
        if res.returncode != 0:
            self.assertIn('usage', res.stderr.strip())
        else:
            self.fail('This test should error')


class DeleteSet(MzidWithDB):
    command = 'deletesets'
    infilename = 'target_pg.tsv'
    dbfn = 'target_psms.sqlite'

    def test_deleteset(self):
        set_to_del = 'Set1'
        sql = """SELECT MIN(rownr), MAX(rownr) FROM psmrows"""
        minrow, maxrow = self.get_values_from_db(self.workdb, sql).fetchone()
        exprownr = 0
        with open(self.infile) as fp:
            head = next(fp).strip('\n').split('\t')
            for line in fp:
                line = {h: l for h, l in zip(head, line.strip('\n').split('\t'))}
                if line['Biological set'] != set_to_del:
                    exprownr += 1
        self.run_command(['--dbfile', self.workdb, '--setnames', 'Set1'])
        newrownr = 0
        with open(self.resultfn) as fp:
            head = next(fp).strip('\n').split('\t')
            for line in fp:
                line = {h: l for h, l in zip(head, line.strip('\n').split('\t'))}
                self.assertNotEqual(line['Biological set'], set_to_del)
                newrownr += 1
        self.assertEqual(newrownr, exprownr)
        newmin, newmax = self.get_values_from_db(self.workdb, sql).fetchone()
        self.assertEqual(minrow, newmin)
        self.assertGreater(maxrow, newmax)
        self.assertEqual(exprownr, newmax + 1)
        sql = """SELECT COUNT(set_name) FROM biosets WHERE set_name='Set1'"""
        self.assertFalse(self.get_values_from_db(self.workdb, sql).fetchone()[0])


class TestIsoSummarize(basetests.MzidTSVBaseTest):
    """Tests producing PSM ratios, not actually summarizing"""
    suffix = '_ratio_isobaric.txt'
    command = 'isosummarize'
    infilename = 'set1_target_pg.tsv'

    def test_mediansweep(self):
        result = self.run_command(['--isobquantcolpattern', 'plex',
            '--mediansweep'])
        self.do_check(0, result.stdout, ratiomethod='sweep')

    def test_keep_zero_NA_psms_denomcolpattern(self):
        result = self.run_command(['--isobquantcolpattern', 'plex',
            '--denompatterns', '_126', '_131', '--keep-psms-na-quant'])
        self.do_check(0, result.stdout, keep_na_quant=True)

    def test_summarize_avg(self):
        result = self.run_command(['--isobquantcolpattern', 'plex',
            '--denompatterns', '_126', '_131', '--summarize-average'])
        self.do_check(0, result.stdout)

    def test_denomcolpattern_regex(self):
        result = self.run_command(['--isobquantcolpattern', 'plex', 
            '--denompatterns', '_1[23][61]'])
        self.do_check(0, result.stdout)

    def get_denominator(self, line, method, denom_ch):
        if method == 'denoms':
            denomvals = [float(line[ch]) for ch in denom_ch if line[ch] != 'NA']
            return sum(denomvals) / len(denomvals)
        elif method == 'sweep':
            return median([float(line[ch]) for ch in line.keys() if line[ch] != 'NA'])

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
                 ratiomethod='denoms', resultch=False, keep_na_quant=False):
        channels = ['tmt10plex_126'] + [x.format('tmt10plex_1', y+27) for x in ['{}{}C', '{}{}N'] for y in range(4)] + ['tmt10plex_131']
        resultch = ['ratio_{}'.format(x) for x in channels]
        denom_ch = [channels[0], channels[-1]]
        if normalize:
            ch_medians = self.check_normalize_medians(channels, denom_ch,
                                                      minint, stdout,
                                                      medianpsms)
        results = [x for x in self.get_values(resultch)]
        resultlinenums = [x[0][0] for x in results]
        for line_num, in_line in enumerate(self.get_infile_lines()):
            in_line.update({ch: in_line[ch]
                            if in_line[ch] != 'NA' and
                            float(in_line[ch]) > minint else 'NA'
                            for ch in channels})
            denom = self.get_denominator({ch: in_line[ch] for ch in channels}, ratiomethod, denom_ch)
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
            if not keep_na_quant and 'NA' in exp_line:
                self.assertNotIn(line_num, resultlinenums)
            else:
                self.assertIn(line_num, resultlinenums)
                nextres = results.pop(0)
                resultline = [x[2] for x in nextres]
                self.assertEqual(resultline, exp_line)




class TestIsoFeatSummarize(basetests.MzidTSVBaseTest):
    suffix = '_ratio_isobaric.txt'
    command = 'isosummarize'
    infilename = 'set1_target_pg.tsv'
    channels = ['tmt10plex_{}'.format(x) for x in ['126', '127N', '127C',
                                                   '128N', '128C', '129N',
                                                   '129C', '130N', '130C',
                                                   '131']]
    nopsms = ['{} - # quanted PSMs'.format(ch) for ch in channels]

    def test_isoquant(self):
        options = ['--featcol', '14', '--isobquantcolpattern', 'tmt10plex',
                   '--denompatterns', '_126']
        self.run_command(options)
        self.isoquant_check(os.path.join(self.fixdir, 'proteins_quantonly.txt'),
            'Master protein(s)', self.channels, self.nopsms)

    def test_normalized_isoquant(self):
        options = ['--featcol', '11', '--isobquantcolpattern', 'tmt10plex',
                   '--denompatterns', '_126', '--keep-psms-na-quant']
        self.run_command(options)
        self.isoquant_check(os.path.join(self.fixdir, 'isosum_charge_column.txt'), # 'proteins_isosum_column.txt'),
            'Charge', self.channels, self.nopsms)
