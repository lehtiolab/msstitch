import os
from numpy import polyfit
from math import log

from tests.integration import basetests


class TestPSM2Peptable(basetests.BaseTest):
    command = 'peptides'
    suffix = '_peptable.tsv'
    infilename = 'set1_target_pg.tsv'
    channels = ['tmt10plex_{}'.format(x) for x in ['126', '127N', '127C',
                                               '128N', '128C', '129N',
                                               '129C', '130N', '130C',
                                               '131']]
    nopsms = ['{} - # quanted PSMs'.format(ch) for ch in channels]

    def check(self, fncol):
        psms = {}
        for psm in self.tsv_generator(self.infile[0]):
            score = str(float(psm['percolator svm-score']))
            try:
                psms[psm['Peptide']][score] = psm
            except KeyError:
                psms[psm['Peptide']] = {score: psm}
        wildcardstrip = 'tmt10plex_'
        striplist = ['#SpecFile', 'Files/scans for peptide', 'MS1 area', 
                'MS1 area (highest of all PSMs)', 'q-value (linear modeled)',
                'Amount fully quanted PSMs']
        for peptide in self.tsv_generator(self.resultfn):
            for newkey, oldkey in zip(['Peptide', 'peptide q-value', 'Protein'],
                                      ['Peptide sequence', 'q-value', 'Protein(s)']):
                peptide[newkey] = peptide.pop(oldkey)
            mapping_psms = psms[peptide['Peptide']]
            toppsm = mapping_psms[str(max([float(score) for score
                                           in mapping_psms.keys()]))]
            testpeptide = {k: v for k, v in peptide.items()
                           if k not in striplist and not k.startswith(wildcardstrip)}
            testpsm = {k: v for k, v in toppsm.items()
                       if k not in striplist and not k.startswith(wildcardstrip)}
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
                             set(['{}_{}'.format(psm[fnfield], psm['SpecID'])
                                  for psm in mapping_psms.values()]))
        self.check_modelqvals()

    def check_modelqvals(self):
        scores, qvalues = [], []
        score, fdr, qthres = 'percolator svm-score', '^q-value', 1e-5
        for line in self.tsv_generator(self.resultfn):
            if float(line[fdr[1:]]) > qthres:
                scores.append(float(line[score]))
                qvalues.append(log(float(line[fdr[1:]]), 10))
        slope, intercept = polyfit(scores, qvalues, deg=1)
        for line in self.tsv_generator(self.resultfn):
            self.assertEqual(float(line[fdr[1:] + ' (linear modeled)']),
                                   10 ** (float(line[score]) *
                                          slope + intercept))

    def check_modelqvals_none_included(self):
        for line in self.tsv_generator(self.resultfn):
            self.assertEqual(line[fdr[1:] + ' (linear modeled)'], 'NA')

    def test_psm2peptable_normalized(self):
        options = ['--spectracol', '1', '--isobquantcolpattern',
                   'tmt10plex', '--scorecolpattern', 'svm',
                   '--ms1quantcolpattern', 'MS1', 
                   '--modelqvals', '--qvalthreshold', '1e-5',
                   '--denompatterns', '126', '--median-normalize',
                   ]
        self.run_command(options)
        self.check(1)
        self.isoquant_check(os.path.join(self.fixdir, 'target_pep_quant_norm.tsv'),
                'Peptide sequence', self.channels, self.nopsms)

    def test_psm2peptable(self):
        fncol = 1
        options = ['--spectracol', str(fncol), '--isobquantcolpattern',
                   'tmt10plex', '--scorecolpattern', 'svm',
                   '--denompatterns', '126',
                   '--ms1quantcolpattern', 'MS1', 
                   '--modelqvals', '--qvalthreshold', '1e-5',
                   '--minpepnr', str(4),
                   ]
        self.run_command(options)
        self.check(fncol)
        self.isoquant_check(os.path.join(self.fixdir, 'target_peptides.tsv'),
                'Peptide sequence', self.channels, self.nopsms)

    def test_average_summarizing(self):
        options = ['--isobquantcolpattern', 'tmt10plex',
                '--scorecolpattern', 'svm',
                   '--denompatterns', '126',
                   '--summarize-average',
                   '--ms1quantcolpattern', 'MS1', 
                   '--modelqvals', '--qvalthreshold', '1e-5',
                   '--minpepnr', str(4),
                   ]
        self.run_command(options)
        self.check(1)
        self.isoquant_check(os.path.join(self.fixdir, 'target_peptides.tsv'),
                'Peptide sequence', self.channels, self.nopsms)

    def test_no_spectracol_mediansweep_keep_na(self):
        options = ['--isobquantcolpattern', 'tmt10plex',
                '--scorecolpattern', 'svm',
                   '--mediansweep',
                   '--keep-psms-na-quant',
                   '--ms1quantcolpattern', 'MS1',
                   '--modelqvals', '--qvalthreshold', '1e-5',
                   '--minpepnr', str(4),
                   ]
        self.run_command(options)
        self.check(1)
        self.isoquant_check(os.path.join(self.fixdir, 'target_peptides_sweep.tsv'),
                'Peptide sequence', self.channels, self.nopsms)

    def test_psm2peptable_totalproteome(self):
        options = ['--spectracol', '1', '--isobquantcolpattern',
                   'tmt10plex', '--scorecolpattern', 'svm',
                   '--ms1quantcolpattern', 'MS1',
                   '--modelqvals', '--qvalthreshold', '1e-5',
                   '--minpepnr', str(4),
                   '--denompatterns', '126',
                   '--totalproteome', os.path.join(self.fixdir, 'proteins.txt'),
                   ]
        self.run_command(options)
        self.check(1)
        self.isoquant_check(os.path.join(self.fixdir, 'target_peptides_totalprotnorm_nolog.txt'),
                'Peptide sequence', self.channels, self.nopsms)

    def test_psm2peptable_totalproteome_logiso(self):
        options = ['--spectracol', '1', '--isobquantcolpattern',
                   'tmt10plex', '--scorecolpattern', 'svm',
                   '--ms1quantcolpattern', 'MS1', 
                   '--modelqvals', '--qvalthreshold', '1e-5',
                   '--minpepnr', str(4),
                   '--denompatterns', '126',
                   '--logisoquant',
                   '--totalproteome', os.path.join(self.fixdir, 'proteins_isonorm_log.txt'),
                   ]
        self.run_command(options)
        self.check(1)
        self.isoquant_check(os.path.join(self.fixdir, 'target_peptides_totalprotnorm.txt'),
                'Peptide sequence', self.channels, self.nopsms)



class TestProteinTable(basetests.ProttableTest):
    command = 'proteins'
    infilename = 'target_peptides.tsv'

    def test_denoms(self):
        self.specialoptions = []
        self.dotest_proteintable('^q-value', 'Master protein(s)', 'Protein ID')
        expectedfn = os.path.join(self.fixdir, 'proteins.txt')
        self.check_lines(expectedfn, self.resultfn)

    def test_sweep(self):
        self.specialoptions = ['--keep-psms-na-quant']
        self.dotest_proteintable('^q-value', 'Master protein(s)', 'Protein ID', summarize_method='sweep')
        expectedfn = os.path.join(self.fixdir, 'proteins_sweep.txt')
        self.check_lines(expectedfn, self.resultfn)

    def test_average_summarizing(self):
        self.specialoptions = ['--summarize-average']
        self.dotest_proteintable('^q-value', 'Master protein(s)', 'Protein ID')
        expectedfn = os.path.join(self.fixdir, 'proteins.txt')
        self.check_lines(expectedfn, self.resultfn)

    def test_isonormalize_log(self):
        self.specialoptions = ['--logisoquant', '--median-normalize', '--minint', '0.1']
        self.dotest_proteintable('^q-value', 'Master protein(s)', 'Protein ID')
        expectedfn = os.path.join(self.fixdir, 'proteins_isonorm_log.txt')
        self.check_lines(expectedfn, self.resultfn)

    def test_isonormalize_nolog_sweep(self):
        self.specialoptions = ['--median-normalize']
        self.dotest_proteintable('^q-value', 'Master protein(s)', 'Protein ID', summarize_method='sweep')
        expectedfn = os.path.join(self.fixdir, 'proteins_isonorm_nolog.txt')
        self.check_lines(expectedfn, self.resultfn)

    def test_no_denom_but_intensity(self):
        self.specialoptions = []
        self.dotest_proteintable('^q-value', 'Master protein(s)', 'Protein ID', summarize_method='intensity')
        expectedfn = os.path.join(self.fixdir, 'proteins_intensities.txt')
        self.check_lines(expectedfn, self.resultfn)

    def test_intensity_normalize(self):
        self.specialoptions = ['--median-normalize']
        res = self.dotest_proteintable('^q-value', 'Master protein(s)', 'Protein ID', summarize_method='intensity', should_error=True)
        if res.returncode != 0:
            self.assertEqual(res.stdout.strip(), 
                    'Cannot do median-centering on intensity values, exiting')
        else:
            self.fail('This test should error due to an invalid combination of '
                    '--median-normalize and --medianintensity')


class TestGenenameTable(basetests.ProttableTest):
    command = 'genes'
    infilename = 'target_peptides.tsv'

    def test(self):
        self.specialoptions = [
                '--targetfasta', os.path.join(self.basefixdir, 'ens99_small.fasta'),
                '--decoyfasta', os.path.join(self.fixdir, 'protrev_ens99_small.fasta')]
        self.dotest_proteintable('^q-value', 'Gene Name', 'Gene Name')
        expectedfn = os.path.join(self.fixdir, 'genenames.txt')
        self.check_lines(expectedfn, self.resultfn)


class TestENSGTable(basetests.ProttableTest):
    command = 'ensg'
    infilename = 'target_peptides.tsv'

    def test(self):
        self.specialoptions = [
                '--targetfasta', os.path.join(self.basefixdir, 'ens99_small.fasta'),
                '--decoyfasta', os.path.join(self.fixdir, 'protrev_ens99_small.fasta')]
        self.dotest_proteintable('svm', 'Gene ID', 'Gene ID')
        expectedfn = os.path.join(self.fixdir, 'ensg.txt')
        self.check_lines(expectedfn, self.resultfn)
