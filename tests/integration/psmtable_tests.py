import os
import re
import sqlite3
import subprocess
from lxml import etree
from Bio import SeqIO
from statistics import median

from app.dataformats import psms as psmh
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
    spectracol = 1
    protkey = 'Protein'
    pepkey = 'Peptide'
    scorekey = 'MSGFScore'
    rt_key = 'Retention time(min)'
    ion_key = 'Ion mobility(Vs/cm2)'
    mc_key = 'missed_cleavage'
    expected_db = 'target_psms.sqlite'
    expected_psms = 'target_pg.tsv'
    timsinfilename ='few_spec_timstof.tsv' 
    timsdbfn = 'spectra_lookup_timstof.sqlite'
    # For not-all-have-isoquant test:
    spec_id_to_rm = '1_controllerType=0 controllerNumber=1 scan=10029'

    """DB and PSM table contain:
    - ENSEMBL proteins
    - a Uniprot swiss protein
    - A self-annotated protein
    - A non-annotated (only peptide) proteins
    """

    def test_build_full_psmtable(self):
        minpif = '0.4'
        fastafn = os.path.join(self.basefixdir, 'ens99_small.fasta')
        options = ['--dbfile', self.workdb,
                '--addbioset', '--genes', '--proteingroup', '--ms1quant', '--isobaric',
                '--fasta', fastafn, '--min-precursor-purity', minpif]
        if self.mc_key:
            options.append('--addmiscleav')
        self.run_command(options)
        self.check_db_fasta(fastafn)
        self.check_addspec_bioset()
        if self.mc_key:
            self.check_miscleav()
        self.check_pg()
        self.check_quanttsv(minpif)
        self.check_addgenes()

    def test_ionmobility_db_inmem(self):
        self.infilename = self.timsinfilename
        self.dbfn = self.timsdbfn
        self.setUp()
        options = ['--dbfile', self.workdb, '--in-memory', '--spectracol', f'{self.spectracol}',
                '--addmiscleav', '--addbioset']
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
        fields = ['Biological set', self.rt_key, 'Ion injection time(ms)', self.ion_key]
        expected_values = self.process_dbvalues_both(self.workdb, sql,
                                                     [1, 2, 3, 4], fields)
        self.check_results_sql(fields, expected_values)
        if self.mc_key:
            self.check_miscleav()

    def test_no_isoquant_in_some_rows(self):
        '''Test when not all rows of PSMs have isoquant data, e.g. when some are
        CID.'''
        # delete isoquant for first scan in PSM table
        db = sqlite3.connect(self.workdb)
        db.execute(f'DELETE FROM isobaric_quant WHERE spectra_id="{self.spec_id_to_rm}"')
        db.execute(f'DELETE FROM precursor_ion_fraction WHERE spectra_id="{self.spec_id_to_rm}"')
        db.commit()
        self.test_build_full_psmtable()

    def check_miscleav(self):
        for val, exp in zip(self.get_values([self.mc_key]), self.get_values([self.pepkey], self.infile[0])):
            exp = re.sub(r'[^A-Za-z]', '', exp[0][0])[:-1]
            self.assertEqual(int(val[0][0]), exp.count('K') + exp.count('R') - exp.count('KP') - exp.count('RP'))

    def check_addspec_bioset(self):
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
        fields = ['Biological set', self.rt_key, 'Ion injection time(ms)', self.ion_key,
                  'Precursor ion fraction']
        expected_values = self.process_dbvalues_both(self.workdb, sql,
                                                     [1, 2, 3, 4, 5], fields)
        self.check_results_sql(fields, expected_values)

    def check_quanttsv(self, minpif):
        sql = ('SELECT pr.rownr, ic.channel_name, '
                'CASE WHEN pif.pif > {} THEN iq.intensity ELSE "NA" END AS intensity '
                'FROM psmrows AS pr JOIN psms USING(psm_id) '
                'LEFT OUTER JOIN precursor_ion_fraction AS pif USING(spectra_id) '
                'JOIN isobaric_quant AS iq USING(spectra_id) '
                'JOIN isobaric_channels AS ic USING(channel_id) '.format(float(minpif)))
        fields = ['tmt10plex_{}'.format(ch) for ch in ['126', '127N', '127C',
                                                       '128N', '128C', '129N',
                                                       '129C', '130N', '130C',
                                                       '131']]
        expected_values = self.rowify(self.get_values_from_db(self.workdb, sql))
        self.check_results_sql(fields, expected_values)
        sql = ('SELECT pr.rownr, pq.intensity '
               'FROM psmrows AS pr JOIN psms USING(psm_id) '
               'LEFT OUTER JOIN ms1_align USING(spectra_id) '
               'LEFT OUTER JOIN ms1_quant AS pq USING(feature_id)')
        expected_values = self.rowify(self.get_values_from_db(self.workdb, sql))
        self.check_results_sql(['MS1 area'], expected_values)

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
                    rsp = rd.split()[1:]
                    nextterm = [i for i,x in enumerate(rsp) if '=' in x][0]
                    desc = ' '.join(rsp[:nextterm])
                elif 'msstitch_fake_gene' in rd:
                    # special case fake fasta record for non-standard gene 
                    six = False
                    desc = rd.split()[-1]
                elif 'msstitch_fake_onlypeptide' in rd:
                    # special fake fasta record for unannotated peptide
                    six, desc = False, 'NA'
                exp_proteins[rec.id] = {
                        'seq': rec.seq,
                        'gene': gene,
                        'desc': desc,
                        'symbol': rd[six: rd.index(' ', six)] if six else 'NA',
                        }
        self.check_db_base(exp_proteins)
        sql = ('SELECT ps.protein_acc, ps.sequence, g.gene_acc, aid.assoc_id, '
                'pd.description '
                'FROM proteins AS p '
                'LEFT OUTER JOIN ensg_proteins USING(pacc_id) '
                'LEFT OUTER JOIN genes AS g USING(gene_id) '
                'LEFT OUTER JOIN genename_proteins AS gnp USING(pacc_id) '
                'LEFT OUTER JOIN associated_ids AS aid USING(gn_id) '
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
                self.assertEqual(exp_proteins[prot]['gene'], gene or 'NA')
                self.assertEqual(exp_proteins[prot]['symbol'], aid or 'NA')
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

    def check_results_sql(self, checkfields, expected_values):
        for resultvals, exp_vals in zip(self.get_values(checkfields),
                                        expected_values):
            for field, resultval, expectval in zip(checkfields, resultvals, exp_vals):
                if field:
                    self.assertEqual([str(x) if x is not None else 'NA'
                                      for x in expectval],
                                     [str(x) for x in resultval])

    def get_expected_psms(self):
        header = self.get_tsvheader(self.infile[0])
        prot_ix = header.index(self.protkey)
        seq_ix = header.index(self.pepkey)
        score_ix = header.index(self.scorekey)
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
        psmhead = psmh.PSMTableHeader([])
        with open(fn) as fp:
            header = next(fp).strip().split('\t')
            master_ix = header.index(psmhead.HEADER_MASTER_PROT)
            pgcontent_ix = header.index(psmhead.HEADER_PG_CONTENT)
            pgamount_ix = header.index(psmhead.HEADER_PG_AMOUNT_PROTEIN_HITS)
            for line in fp:
                line = line.strip().split('\t')
                yield {'master': line[master_ix],
                       'content': line[pgcontent_ix],
                       'amount': line[pgamount_ix],
                       }

    def check_pglup(self, sql, keyfun, valfun):
        result = {keyfun(x): valfun(x) for x in
                  self.get_values_from_db(self.workdb, sql)}
        exp_file = os.path.join(self.fixdir, self.expected_db)
        expected = {keyfun(x): valfun(x) for x in
                    self.get_values_from_db(exp_file, sql)}
        for key, value in result.items():
            self.assertIn(key, expected.keys())
            if type(value) == tuple:
                for v, e in zip(value, expected[key]):
                    try:
                        self.assertEqual(v, e)
                    except AssertionError:
                        # I dont understand why but sometimes the expected value
                        # is a "long float" 0.12345000000001, and sometimes it's 0.12345
                        # making this test flaky when comparing protein_score
                        # which is a REAL (float) for sage
                        self.assertAlmostEqual(v, e)
            else:
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
        expected = self.parse_proteingroups(os.path.join(self.fixdir, self.expected_psms))
        for res, exp in zip(result, expected):
            self.assertEqual(set(res['master'].split(';')),
                             set(exp['master'].split(';')))
            self.assertEqual(res['amount'], exp['amount'])
            rescontent = res['content'].split(';')
            expcontent = exp['content'].split(';')
            self.assertEqual(set(rescontent), set(expcontent))


    def check_addgenes(self):
        for line in self.get_values(['Gene ID', 'Gene Name', 'Description', self.protkey]):
            genes = line[0][1].split(';')
            assoc_ids = line[1][1].split(';')
            descriptions = ['{}]'.format(x).replace(']]', ']')
                            for x in line[2][1].split('];')]
            proteins = [x.split('(')[0] for x in line[3][1].split(';')]
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


class TestPSMTableSage(TestPSMTable):
    infilename = 'target.sage.tsv'
    spectracol = 5
    protkey = 'proteins'
    pepkey = 'peptide'
    scorekey = 'sage_discriminant_score'
    rt_key = False
    ion_key = False
    mc_key = False
    expected_db = 'target_psms.sage.sqlite'
    expected_psms = 'target_pg.sage.tsv'
    timsinfilename ='few_spec_timstof.sage.tsv' 
    spec_id_to_rm = 'controllerType=0 controllerNumber=1 scan=10192'


class TestPercoTSV(basetests.MzidTSVBaseTest):
    command = 'perco2psm'
    suffix = '_fdr.tsv'
    infilename = 'few_spectra.tsv'
    mzidfn = 'few_spectra.mzid'
    percofn = 'perco.xml'
    expected_fn = 'few_spectra.tsv_fdr.tsv'
    qval_psms = 'qvality_psms.txt'
    qval_peps = 'qvality_peps.txt'


    def test_conffilt(self):
        threshold = 0.3
        options = ['--perco', os.path.join(self.fixdir, self.percofn), '--filtpep', 
            str(threshold), '--filtpsm', str(threshold)]
        checkfields = ['percolator svm-score', 'PSM q-value', 'peptide q-value',
            'PSM PEP', 'peptide PEP', 'TD']
        if self.mzidfn:
            options.extend(['--mzids', os.path.join(self.fixdir, self.mzidfn)])
            checkfields.append('ExplainedIonCurrentRatio')
        self.run_command(options)
        with open(os.path.join(self.fixdir, self.expected_fn)) as fp:
            header = next(fp).strip().split('\t')
            expected = [line.strip().split('\t') for line in fp]
        expected = [{field: line[i] for i, field in enumerate(header)} for line in expected]
        expected = [line for line in expected if float(line['PSM q-value']) < threshold and
            float(line['peptide q-value']) < threshold]
        # check how many PSMs made the cut:
        self.assertEqual(len(expected),  len([x for x in self.get_values(checkfields)]))
        for res, exp in zip(self.get_values(checkfields), expected):
            for i, field in enumerate(checkfields):
                self.assertEqual(field, res[i][0])
                self.assertEqual(exp[field], res[i][1])

    def test_add_fdr(self):
        options = ['--perco', os.path.join(self.fixdir, self.percofn)]
        checkfields = ['percolator svm-score', 'PSM q-value', 'peptide q-value',
            'PSM PEP', 'peptide PEP', 'TD']
        if self.mzidfn:
            options.extend(['--mzids', os.path.join(self.fixdir, self.mzidfn)])
            checkfields.append('ExplainedIonCurrentRatio')
        self.run_command(options)
        with open(os.path.join(self.fixdir, self.expected_fn)) as fp:
            header = next(fp).strip().split('\t')
            expected = [line.strip().split('\t') for line in fp]
        expected = [{field: line[i] for i, field in enumerate(header)} for line in expected]
        self.assertEqual(len(expected),  len([x for x in self.get_values(checkfields)]))
        for res, exp in zip(self.get_values(checkfields), expected):
            for i, field in enumerate(checkfields):
                self.assertEqual(field, res[i][0])
                self.assertEqual(exp[field], res[i][1])

    def test_use_qvality_tsv(self):
        '''Get FDR data from recalculation in qvality. Qvality input is prepared
        as a TSV with header, for target and decoy (so qvality -d) . In this case 
        I have simply taken the percolator svm/qval/pep scores from perco.xml and 
        added 0.1 to each PSM (not peptide) q-value using bash, PSMs (double backslashes
        for python):
        cat <(echo $'Score\tPEP\tq-value') \
                <(grep -A3 '<psm p' perco.xml| grep -v psm| grep -v '\\-\\-' | \
                sed -E 's/\\s*<[\\/]*[a-z_]+>//g' | paste - - - | \
                gawk -v FS='\t' -v OFS='\t' '{print $1, $3, $2+0.1}') > qvality_psms.txt 

        Peptides:
        cat ...as above \
                <(grep -A3 '<peptide p' perco.xml| grep -v peptide| grep -v '\\-\\-' | \
                sed -E 's/\\s*<[\\/]*[a-z_]+>//g' | paste - - - | \
                gawk -v FS='\t' -v OFS='\t' '{print $1, $3, $2}') > qvality_peps.txt 

        '''
        options = ['--perco', os.path.join(self.fixdir, self.percofn),
                '--qvalitypsms', os.path.join(self.fixdir, self.qval_psms),
                '--qvalitypeps', os.path.join(self.fixdir, self.qval_peps)]
        if self.mzidfn:
            options.extend(['--mzids', os.path.join(self.fixdir, self.mzidfn)])
        self.run_command(options)
        checkfields = ['percolator svm-score', 'PSM q-value', 'peptide q-value',
            'PSM PEP', 'peptide PEP', 'TD']
        with open(os.path.join(self.fixdir, self.expected_fn)) as fp:
            header = next(fp).strip().split('\t')
            expected = [line.strip().split('\t') for line in fp]
        expected = [{field: line[i] for i, field in enumerate(header)} for line in expected]
        self.assertEqual(len(expected),  len([x for x in self.get_values(checkfields)]))
        for res, exp in zip(self.get_values(checkfields), expected):
            for i, field in enumerate(checkfields):
                if field == 'PSM q-value':
                    self.assertAlmostEqual(float(exp[field]) + 0.1, float(res[i][1]), places=6)
                else:
                    self.assertEqual(exp[field], res[i][1])
                self.assertEqual(field, res[i][0])


class TestPercoTSVTIMS(basetests.MzidTSVBaseTest):
    command = 'perco2psm'
    suffix = '_fdr.tsv'
    infilename = 'few_spec_timstof.tsv'
    dbfn = 'spectra_lookup_timstof.sqlite'

    def test_add_fdr_timstof(self):
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
                self.assertEqual(field, res[i][0])
                self.assertEqual(exp[field], res[i][1])


class TestPercoTSVSage(TestPercoTSV):
    infilename = 'few_spectra.sage.tsv'
    mzidfn = False
    percofn = 'perco.sage.xml'
    expected_fn = 'few_spectra.sage.fdr.tsv'
    qval_psms = 'qvality_psms.sage.txt'
    qval_peps = 'qvality_peps.sage.txt'


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
        options = ['--splitcol', '31']
        self.run_command(options)
        resultfiles = [os.path.join(self.workdir, '{}.tsv'.format(setname))
                       for setname in setnames]
        for resultfn in resultfiles:
            for line in self.get_all_lines(resultfn):
                self.assertIn(line, self.expectlines)


class TestSeqMatchFastaDB(basetests.MzidTSVBaseTest):
    command = 'seqmatch'
    infilename = 'few_spectra.tsv'
    suffix = '_seqmatch.txt'

    def check_peps_in_out(self, options, seqs, matching=True):
        matchcol = 'TESTMATCH'
        self.run_command([*options, '--dbfile', 'seqs.db', '--matchcolname', matchcol])
        pepfound = False
        seqs_nrmap = {x: ix for ix,x in enumerate(seqs)}
        with open(self.resultfn) as fp:
            header = next(fp).strip().split('\t')
            out_linecount = 0
            for line in fp:
                out_linecount += 1
                psm = {k: v for k,v in zip(header, line.strip().split('\t'))}
                seq = re.sub(r'[^A-Z]', '', psm[self.pepkey])
                if seq in seqs_nrmap and matching:
                    self.assertEqual(psm[matchcol], f'acc_{seqs_nrmap[seq]}')
                else:
                    self.assertEqual(psm[matchcol], 'No match')
        with open(self.infile[0]) as fp:
            header = next(fp).strip().split('\t')
            in_linecount = 0
            for line in fp:
                in_linecount += 1
        self.assertEqual(in_linecount, out_linecount)

    def test_noflags(self):
        seqs = ['DGTDVLR']
        basetests.create_db(seqs, mapaccession=True)
        options = ['--dbfile', 'seqs.db']
        self.check_peps_in_out(options, seqs)

    def test_ntermwildcards_inmem(self):
        seqs = ['XXDGTDVLR']
        seqs_to_filter = ['DGTDVLR']
        basetests.create_db(seqs, reverse=True, mapaccession=True)

        # Find peptide
        max_falloff = 2
        options = ['--insourcefrag', str(max_falloff), '--in-memory']
        self.check_peps_in_out(options, seqs_to_filter)

# FIXME
#        # Now failing, too short falloff
#        max_falloff = 1
#        options = ['--insourcefrag', str(max_falloff)]
#        self.check_peps_in_out(options, seqs_to_filter)

    def test_deamidate(self):
        # deamidation: N -> D
        seqs = ['NGTDVIR']
        seqs_to_filter = ['DGTDVLR']
        basetests.create_db(seqs, mapaccession=True)
        options = ['--deamidate']
        self.check_peps_in_out(options, seqs_to_filter)

    def test_fullprotein(self):
        '''Testing if PSMs with sequence match (except Isoleucine) are removed'''
        seqs = ['XXXXXXXXXXDGTDVIRYYYYYYY']
        seqs_to_filter = ['DGTDVLR']
        basetests.create_db(seqs, fullprotein=True, minlen=7, mapaccession=True)
        # Find peptide
        options = ['--fullprotein', '--minlen', '7']
        self.check_peps_in_out(options, seqs_to_filter)

    def test_fullprotein_minlenwrong(self):
        '''Do not filter since even with matching 9-aa peptide, the 7th of the PSM does not match'''
        seqs = ['XXXXXXXXXXDGTDVLYYYYYYY']
        seqs_to_filter = ['DGTDVLR']
        basetests.create_db(seqs, fullprotein=True, minlen=7, mapaccession=True)
        # Find peptide
        options = ['--fullprotein', '--minlen', '7']
        self.check_peps_in_out(options, seqs_to_filter, matching=False)


class TestSeqMatchFastaDBSage(TestSeqMatchFastaDB):
    infilename = 'few_spectra.sage.tsv'
    pepkey = 'peptide'


class TestRemoveDupPSM(basetests.MzidTSVBaseTest):
    command = 'deduppsms'
    infilename = 'few_spectra_duplicate.tsv'
    suffix = '_dedup.txt'

    def get_line_count(self, fn):
        with open(fn) as fp:
            return sum(1 for _line in fp)

    def test_nodup_noremove(self):
        options = []
        self.infile = [os.path.join(self.fixdir, 'few_spectra.tsv')]
        self.run_command(options)
        self.assertEqual(self.get_line_count(self.infile[0]), self.get_line_count(self.resultfn))

    def test_removedup(self):
        '''There are two duplicated PSMs in the input file, which is derived from
        the few_spectra.tsv file'''
        options = []
        self.run_command(options)
        self.assertEqual(self.get_line_count(self.infile[0]) - 2, self.get_line_count(self.resultfn))
        checkfile = os.path.join(self.fixdir, 'few_spectra.tsv')
        self.assertEqual(os.path.getsize(checkfile), os.path.getsize(self.resultfn))

    def test_removedup_seqcol(self):
        '''Some PSMs are already duplicate by Scan Number since they have 
        two equal scoring solutions from search engine by sequence'''
        self.infile = [os.path.join(self.fixdir, 'few_spectra.tsv')]
        options = ['--peptidecolpattern', 'ScanNum']
        self.run_command(options)
        self.assertEqual(self.get_line_count(self.infile[0]) - 3, self.get_line_count(self.resultfn))


class TestSeqFilt(basetests.MzidTSVBaseTest):
    command = 'seqfilt'
    infilename = 'few_spectra.tsv'
    suffix = '_filtseq.txt'
    specidkey = 'SpecID'

    def test_noflags_inmem(self):
        seqs = ['DGTDVLR']
        basetests.create_db(seqs)
        options = ['--dbfile', 'seqs.db', '--in-memory']
        self.run_command(options)
        self.check_peps_in_out(seqs, matching=True)

    def check_peps_in_out(self, seqs, matching):
        '''Check if peptides in seqs are contained in new PSM tables.
        matching is 0/1 - True: peptides passed match, should NOT be found anymore
                          False: peptides do not match, should be found in DB
        '''
        pepfound = False
        with open(self.resultfn) as fp:
            header = next(fp).strip().split('\t')
            filtpsms = {}
            for line in fp:
                psm = {k: v for k,v in zip(header, line.strip().split('\t'))}
                seq = re.sub(r'[^A-Z]', '', psm[self.pepkey])
                scan = psm[self.specidkey].split('=')[-1]
                # Use scan_seq since sometimes have two solutions to same scan at equal score
                filtpsms[f'{scan}_{seq}'] = seq
                if seq in seqs:
                    pepfound = True
                    break
        if matching:
            self.assertFalse(pepfound)
        else:
            self.assertTrue(pepfound)
        with open(self.infile[0]) as fp:
            allpsms = {}
            header = next(fp).strip().split('\t')
            for line in fp:
                psm = {k: v for k,v in zip(header, line.strip().split('\t'))}
                seq = re.sub(r'[^A-Z]', '', psm[self.pepkey])
                scan = psm[self.specidkey].split('=')[-1]
                allpsms[f'{scan}_{seq}'] = seq
                if seq in seqs:
                    pepfound = True
                    break
                else:
                    self.assertEqual(seq, filtpsms[f'{scan}_{seq}'])
        self.assertTrue(pepfound)

    def test_ntermwildcards(self):
        seqs = ['XXDGTDVLR']
        seqs_to_filter = ['DGTDVLR']
        basetests.create_db(seqs, reverse=True)

        # Find peptide
        max_falloff = 2
        options = ['--dbfile', 'seqs.db', '--insourcefrag', str(max_falloff)]
        self.run_command(options)
        self.check_peps_in_out(seqs_to_filter, matching=True)

        # Now failing, too short falloff
        max_falloff = 1
        options = ['--dbfile', 'seqs.db', '--insourcefrag', str(max_falloff)]
        self.run_command(options)
        self.check_peps_in_out(seqs_to_filter, matching=False)

    def test_deamidate(self):
        # deamidation: N -> D
        seqs = ['NGTDVIR']
        seqs_to_filter = ['DGTDVLR']
        basetests.create_db(seqs)
        options = ['--dbfile', 'seqs.db', '--deamidate']
        self.run_command(options)
        self.check_peps_in_out(seqs_to_filter, matching=True)

    def test_fullprotein(self):
        '''Testing if PSMs with sequence match (except Isoleucine) are removed'''
        seqs = ['XXXXXXXXXXDGTDVIRYYYYYYY']
        seqs_to_filter = ['DGTDVLR']
        basetests.create_db(seqs, fullprotein=True, minlen=6)
        # Find peptide
        options = ['--dbfile', 'seqs.db', '--fullprotein', '--minlen', '6']
        self.run_command(options)
        self.check_peps_in_out(seqs_to_filter, matching=True)

    def test_fullprotein_minlenwrong(self):
        '''Do not filter since even with matching 6-aa peptide, the 7th of the PSM does not match'''
        seqs = ['XXXXXXXXXXDGTDVLYYYYYYY']
        seqs_to_filter = ['DGTDVLR']
        basetests.create_db(seqs, fullprotein=True, minlen=6)
        # Find peptide
        options = ['--dbfile', 'seqs.db', '--fullprotein', '--minlen', '9']
        self.run_command(options)
        self.check_peps_in_out(seqs_to_filter, matching=False)


class TestSeqFiltSage(TestSeqFilt):
    infilename = 'few_spectra.sage.tsv'
    pepkey = 'peptide'
    specidkey = 'scannr'


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


class TestDeleteSet(MzidWithDB):
    command = 'deletesets'
    infilename = 'target_pg.tsv'
    dbfn = 'target_psms.sqlite'

    def do_deleteset(self, options):
        set_to_del = 'Set1'
        db = sqlite3.connect(self.workdb)
        # Add peptide which is to be deleted since it will not belong to a remaining
        # set (peptides in set1/2 in DB are identical)
        db.execute('INSERT INTO peptide_sequences(sequence) VALUES("fake_seq")')
        db.commit()
        
        rowsql = 'SELECT MIN(rownr), MAX(rownr) FROM psmrows'
        minrow, maxrow = db.execute(rowsql).fetchone()
        pepsql = 'SELECT COUNT(*) FROM peptide_sequences'
        pep_precount = db.execute(pepsql).fetchone()[0]
        exprownr = 0
        with open(self.infile) as fp:
            head = next(fp).strip('\n').split('\t')
            for line in fp:
                line = {h: l for h, l in zip(head, line.strip('\n').split('\t'))}
                if line['Biological set'] != set_to_del:
                    exprownr += 1
        self.run_command(['--dbfile', self.workdb, '--setnames', set_to_del, *options])
        newrownr = 0
        with open(self.resultfn) as fp:
            head = next(fp).strip('\n').split('\t')
            for line in fp:
                line = {h: l for h, l in zip(head, line.strip('\n').split('\t'))}
                self.assertNotEqual(line['Biological set'], set_to_del)
                newrownr += 1
        self.assertEqual(newrownr, exprownr)
        newmin, newmax = db.execute(rowsql).fetchone()
        self.assertEqual(minrow, newmin)
        self.assertGreater(maxrow, newmax)
        self.assertEqual(exprownr, newmax + 1)
        pep_newcount = db.execute(pepsql).fetchone()[0]
        self.assertGreater(pep_precount, pep_newcount)
        biosql = """SELECT COUNT(set_name) FROM biosets WHERE set_name='Set1'"""
        self.assertFalse(db.execute(biosql).fetchone()[0])

    def test_db_fn(self):
        self.do_deleteset([])

    def test_db_inmem(self):
        self.do_deleteset(['--in-memory'])


class TestDeleteSetSage(TestDeleteSet):
    infilename = 'target_pg.sage.tsv'
    dbfn = 'target_psms.sage.sqlite'


class TestIsoSummarize(basetests.MzidTSVBaseTest):
    """Tests producing PSM ratios, not actually summarizing"""
    suffix = '_ratio_isobaric.txt'
    command = 'isosummarize'
    infilename = 'set1_target_pg.tsv'

    def test_mediansweep(self):
        result = self.run_command(['--isobquantcolpattern', 'plex',
            '--mediansweep'])
        self.do_check(0, result.stdout, ratiomethod='sweep')

    def test_sage(self):
        self.infile = os.path.join(self.fixdir, 'set1_target_pg.sage.tsv')
        result = self.run_command(['--isobquantcolpattern', 'plex',
            '--mediansweep'])
        self.do_check(0, result.stdout, ratiomethod='sweep', search='sage')

    def test_summarize_avg(self):
        denompats = ['_126']
        result = self.run_command(['--isobquantcolpattern', 'plex',
            '--denompatterns', *denompats, '--summarize-average'])
        self.do_check(0, result.stdout, denompats=denompats)

    def test_denomcolpattern_regex(self):
        denompats = ['_1[23][61]']
        result = self.run_command(['--isobquantcolpattern', 'plex', 
            '--denompatterns', *denompats])
        self.do_check(0, result.stdout, denompats=denompats)

    def get_denominator(self, line, method, denom_ch):
        if method == 'denoms':
            denomvals = [float(line[ch]) for ch in denom_ch if line[ch] != 'NA']
            if denomvals == []:
                return 0
            else:
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

    def do_check(self, minint, stdout, normalize=False, medianpsms=None, ratiomethod='denoms',
            resultch=False, denompats=False, search=False):
        channels = ['tmt10plex_126'] + [x.format('tmt10plex_1', y+27) for x in ['{}{}C', '{}{}N'] for y in range(4)] + ['tmt10plex_131']
        if search == 'sage':
            spec_id_fields = ['filename', 'scannr']
        else:
            spec_id_fields = ['#SpecFile', 'SpecID'] 
            
        resultch = spec_id_fields + [f'ratio_{x}' for x in channels]
        denom_ch = []
        if denompats:
            for denompat in denompats:
                denom_ch.extend([x for x in channels if re.search(denompat, x)])
        if normalize:
            ch_medians = self.check_normalize_medians(channels, denom_ch,
                                                      minint, stdout,
                                                      medianpsms)
        results = [x for x in self.get_values(resultch)]
        resultspecids = [f'{x[0][1]}_{x[1][1]}' for x in results]
        for line_num, in_line in enumerate(self.get_infile_lines()):
            in_line.update({ch: in_line[ch]
                            if in_line[ch] != 'NA' and
                            float(in_line[ch]) > minint else 'NA'
                            for ch in channels})
            specid = f'{in_line[spec_id_fields[0]]}_{in_line[spec_id_fields[1]]}'
            denom = self.get_denominator({ch: in_line[ch] for ch in channels},
                    ratiomethod, denom_ch)
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
            self.assertIn(specid, resultspecids)
            nextres = results.pop(0)
            resultline = [x[1] for x in nextres]
            self.assertEqual(resultline[2:], exp_line)


class TestIsoFeatSummarize(basetests.MzidTSVBaseTest):
    suffix = '_ratio_isobaric.txt'
    command = 'isosummarize'
    infilename = 'set1_target_pg.tsv'
    channels = ['tmt10plex_{}'.format(x) for x in ['126', '127N', '127C',
                                                   '128N', '128C', '129N',
                                                   '129C', '130N', '130C',
                                                   '131']]
    nopsms = ['{} - Quanted PSM count'.format(ch) for ch in channels]

    def test_isoquant(self):
        options = ['--featcol', '14', '--isobquantcolpattern', 'tmt10plex',
                   '--denompatterns', '_126']
        self.run_command(options)
        self.isoquant_check(os.path.join(self.fixdir, 'proteins_quantonly.txt'),
            'Master protein(s)', self.channels, self.nopsms)

    def test_sage(self):
        self.infile = os.path.join(self.fixdir, 'set1_target_pg.sage.tsv')
        options = ['--featcol', '5', '--isobquantcolpattern', 'tmt10plex',
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

    def test_split_multi_mapping(self):
        options = ['--featcol', '13', '--isobquantcolpattern', 'tmt10plex',
                   '--denompatterns', '_126', '--split-multi-entries']
        self.run_command(options)
        self.isoquant_check(os.path.join(self.fixdir, 'proteins_quantonly_splitmulti.txt'),
            'Protein', self.channels, self.nopsms)
