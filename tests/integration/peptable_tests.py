import os
from numpy import polyfit
from math import log

from tests.integration import basetests


class TestPSM2Peptable(basetests.PeptableTest):
    command = 'psm2pep'
    suffix = '_peptable.tsv'
    infilename = 'mzidtsv.txt'

    def check(self, fncol):
        psms = {}
        for psm in self.tsv_generator(self.infile[0]):
            score = str(float(psm['percolator svm-score']))
            try:
                psms[psm['Peptide']][score] = psm
            except KeyError:
                psms[psm['Peptide']] = {score: psm}
        striplist = ['fake_ch{}'.format(x) for x in range(8)]
        striplist.extend(['Files/scans for peptide', 'MS1 area',
                          'MS1 area (highest of all PSMs)'])
        for peptide in self.tsv_generator(self.resultfn):
            for newkey, oldkey in zip(['Peptide', 'peptide PEP',
                                       'peptide q-value', 'Protein'],
                                      ['Peptide sequence', 'PEP', 'q-value',
                                       'Protein(s)']):
                peptide[newkey] = peptide.pop(oldkey)
            mapping_psms = psms[peptide['Peptide']]
            toppsm = mapping_psms[str(max([float(score) for score
                                           in mapping_psms.keys()]))]
            testpeptide = {k: v for k, v in peptide.items()
                           if k not in striplist}
            testpsm = {k: v for k, v in toppsm.items()
                       if k not in striplist}
            self.assertEqual(testpeptide, testpsm)
            ms1s = [self.get_float_or_na(psm['MS1 area'])
                    for psm in mapping_psms.values()]
            try:
                topms1 = max(ms1s)
            except TypeError:
                # cannot do max on mix of str and float
                topms1 = max([0 if x == 'NA' else x for x in ms1s])
            self.assertEqual(str(topms1),
                             peptide['MS1 area (highest of all PSMs)'])
            psmfnscans = set(peptide['Files/scans for peptide'].split('; '))
            with open(self.infile[0]) as fp:
                fnfield = next(fp).strip().split('\t')[fncol - 1]
            self.assertEqual(psmfnscans,
                             set(['{}_{}'.format(psm[fnfield], psm['ScanNum'])
                                  for psm in mapping_psms.values()]))

    def test_psm2peptable(self):
        fncol = 1
        options = ['--spectracol', str(fncol), '--isobquantcolpattern',
                   'fake_ch', '--scorecolpattern', 'svm',
                   '--ms1quantcolpattern', 'MS1']
        self.run_command(options)
        self.check(fncol)

    def test_no_spectracol(self):
        options = ['--isobquantcolpattern', 'fake_ch',
                   '--scorecolpattern', 'svm', '--ms1quantcolpattern', 'MS1']
        self.run_command(options)
        self.check(1)


class TestIsoquant(basetests.PeptableTest):
    infilename = 'peptable_small.txt'
    command = 'isoquant'
    suffix = '_isoq.txt'
    channels = ['fake_ch{}'.format(x) for x in range(0, 8)]
    nopsms = ['{} - # quanted PSMs'.format(ch) for ch in channels]

    def test_isoquant(self):
        isotable = os.path.join(self.fixdir, 'peptable_isoquant.txt')
        options = ['--quantfile', isotable, '--isobquantcolpattern', 'fake_ch',
                   '--qaccpattern', 'accession']
        self.run_command(options)
        self.isoquant_check(isotable, 'Peptide sequence', self.channels,
                            self.nopsms)


class TestModelQvals(basetests.PeptableTest):
    command = 'modelqvals'
    infilename = 'peptable_small.txt'
    suffix = '_qmodel.txt'

    def test_modelqvals(self):
        score, fdr = 'percolator svm-score', '^q-value'
        options = ['--scorecolpattern', score, '--fdrcolpattern', fdr]
        self.run_command(options)
        scores, qvalues = [], []
        for line in self.tsv_generator(self.infile[0]):
            if float(line[fdr[1:]]) > 10e-4:
                scores.append(float(line[score]))
                qvalues.append(log(float(line[fdr[1:]]), 10))
        slope, intercept = polyfit(scores, qvalues, deg=1)
        for line in self.tsv_generator(self.resultfn):
            self.assertAlmostEqual(float(line[fdr[1:] + ' (linear modeled)']),
                                   10 ** (float(line[score]) *
                                          slope + intercept))


class TestBuild(basetests.PeptableTest):
    command = 'build'
    infilename = 'built_peptide_table.txt'
    suffix = ''

    def setUp(self):
        super().setUp()
        self.dbfile = os.path.join(self.fixdir, 'peptable_db.sqlite')

    def get_std_options(self):
        cmd = [self.executable, self.command, '-d', self.workdir]
        return cmd

    def check(self, genecentric=False, noncentric=False):
        self.check_values()
        if genecentric:
            self.check_genes()
        elif noncentric:
            pass
        else:
            self.check_proteingroups()
        self.check_isobaric()

    def check_values(self):
        sql = ('SELECT ps.sequence, bs.set_name, '
               'ppq.quant, pf.fdr, pp.pep '
               'FROM peptide_sequences AS ps '
               'JOIN biosets AS bs '
               'JOIN peptide_precur_quanted AS ppq USING(pep_id) '
               'JOIN peptide_fdr AS pf USING(pep_id) '
               'JOIN peptide_pep AS pp USING(pep_id) '
               )
        self.check_build_values(sql, ['MS1 area (highest of all PSMs)',
                                      'q-value', 'PEP'], 'Peptide sequence')

    def check_genes(self):
        # FIXME get a db WITHOUT protein groups for genecentric peptides
        sql = ('SELECT ps.sequence, p.psm_id, "NA", pd.description, '
               'g.gene_acc, aid.assoc_id, "NA" '
               'FROM peptide_sequences AS ps '
               'JOIN psms AS p USING(pep_id) '
               'JOIN protein_psm USING(psm_id) '
               'JOIN prot_desc AS pd USING(protein_acc) '
               'JOIN genes AS g USING(protein_acc) '
               'JOIN associated_ids AS aid USING(protein_acc) '
               'JOIN protein_coverage AS pc USING(protein_acc) '
               )
        self.check_peptide_relations(sql)

    def check_proteingroups(self):
        # FIXME have to also check content proteins, etc
        sql = ('SELECT ps.sequence, p.psm_id, pm.protein_acc, pd.description, '
               'g.gene_acc, aid.assoc_id, pc.coverage '
               'FROM peptide_sequences AS ps '
               'JOIN psms AS p USING(pep_id) '
               'JOIN psm_protein_groups USING(psm_id) '
               'JOIN protein_group_master AS pm USING(master_id) '
               'JOIN prot_desc AS pd USING(protein_acc) '
               'JOIN genes AS g USING(protein_acc) '
               'JOIN associated_ids AS aid USING(protein_acc) '
               'JOIN protein_coverage AS pc USING(protein_acc) '
               )
        self.check_peptide_relations(sql)

    def check_peptide_relations(self, sql):
        expected, psm_id, pep = {}, None, None
        for rec in self.get_values_from_db(self.dbfile, sql):
            try:
                expected[rec[0]]['psms'].add(rec[1])
            except KeyError:
                expected[rec[0]] = {'psms': set([rec[1]]),
                                    'pgroups': set([rec[2]]),
                                    'descriptions': set([rec[3]]),
                                    'genes': set([rec[4]]),
                                    'assoc': set([rec[5]]),
                                    'cover': set([str(rec[6])]),
                                    }
            else:
                expected[rec[0]]['pgroups'].add(rec[2])
                expected[rec[0]]['descriptions'].add(rec[3])
                expected[rec[0]]['genes'].add(rec[4])
                expected[rec[0]]['assoc'].add(rec[5])
                expected[rec[0]]['cover'].add(str(rec[6]))
        for line in self.tsv_generator(self.resultfn):
            self.assertEqual(set(line['Protein(s)'].split(';')),
                             expected[line['Peptide sequence']]['pgroups'])
            self.assertEqual(
                set(line['Description(s)'].split(';')),
                set([y for x in
                     expected[line['Peptide sequence']]['descriptions']
                     for y in x.split(';')]))
            self.assertEqual(set(line['Gene(s)'].split(';')),
                             expected[line['Peptide sequence']]['genes'])
            self.assertEqual(set(line['Associated gene ID(s)'].split(';')),
                             expected[line['Peptide sequence']]['assoc'])
            self.assertEqual(set([line['Coverage(s)']]),
                             expected[line['Peptide sequence']]['cover'])

    def check_isobaric(self):
        sql = ('SELECT ps.sequence, bs.set_name, ch.channel_name, '
               'iq.quantvalue, iq.amount_psms '
               'FROM peptide_sequences AS ps '
               'JOIN biosets AS bs '
               'JOIN peptide_iso_quanted AS iq USING(pep_id) '
               'JOIN pepquant_channels AS ch USING(channel_id) '
               )
        self.check_built_isobaric(sql, 'Peptide sequence')

    def test_proteincentric(self):
        options = ['--fdr', '--pep', '--isobaric', '--precursor',
                   '--dbfile', self.dbfile]
        self.run_command(options)
        self.check()

    def test_genecentric(self):
        options = ['--fdr', '--pep', '--isobaric', '--precursor',
                   '--dbfile', self.dbfile, '--genecentric']
        self.run_command(options)
        self.check(genecentric=True)

    def test_noncentric(self):
        options = ['--fdr', '--pep', '--isobaric', '--precursor',
                   '--dbfile', self.dbfile, '--noncentric']
        self.run_command(options)
        self.check(noncentric=True)
