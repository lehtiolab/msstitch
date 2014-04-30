import re
import os
import basetest
import sqlite3
import yaml


class TestSplitTD(basetest.BaseTestPycolator):
    command = 'splittd'
    infilename = 'percolator_out.xml'
    suffix = ''

    def test_split(self):
        """Tests that splitted files contain equal amount of PSMS
        when compared with expected output, and checks that each psm/peptide
        has correct 'decoy' attribute."""
        self.target = os.path.join(self.workdir,
                                   self.infilename + '_target.xml')
        self.decoy = os.path.join(self.workdir,
                                  self.infilename + '_decoy.xml')
        self.run_pycolator()
        target_expected = os.path.join(self.fixdir, 'splittd_target_out.xml')
        decoy_expected = os.path.join(self.fixdir, 'splittd_decoy_out.xml')
        target_exp_contents = self.read_percolator_out(target_expected)
        decoy_exp_contents = self.read_percolator_out(decoy_expected)
        target_contents = self.read_percolator_out(self.target)
        decoy_contents = self.read_percolator_out(self.decoy)

        self.assertEqual(len(target_contents['psms']),
                         len(target_exp_contents['psms']))
        self.assertEqual(len(target_contents['peptides']),
                         len(target_exp_contents['peptides']))
        self.assertEqual(len(decoy_contents['psms']),
                         len(decoy_exp_contents['psms']))
        self.assertEqual(len(decoy_contents['peptides']),
                         len(decoy_exp_contents['peptides']))
        for feat in ['psms', 'peptides']:
            for el in target_contents[feat]:
                self.assertEqual(
                    el.attrib['{%s}decoy' % target_contents['ns']], 'false')
            for el in decoy_contents[feat]:
                self.assertEqual(
                    el.attrib['{%s}decoy' % decoy_contents['ns']], 'true')


class TestMerge(basetest.BaseTestPycolator):
    command = 'merge'
    infilename = 'splittd_target_out.xml'
    suffix = '_merged.xml'

    def test_merge(self):
        self.multifiles = [os.path.join(self.fixdir, 'splittd_decoy_out.xml')]
        options = ['--multifiles']
        options.extend(self.multifiles)
        self.run_pycolator(options)
        expected = self.read_percolator_out(os.path.join(self.fixdir,
                                                         'percolator_out.xml'))
        result = self.read_percolator_out(self.resultfn)
        self.assertEqual(len(result['psms']), len(result['peptides']))
        self.assertCountEqual(self.get_element_ids(expected['psms'], 'psm_id',
                                                   expected['ns']),
                              self.get_element_ids(result['psms'], 'psm_id',
                                                   result['ns']))
        self.assertCountEqual(self.get_element_ids(expected['peptides'],
                                                   'peptide_id',
                                                   expected['ns']),
                              self.get_element_ids(result['peptides'],
                                                   'peptide_id', result['ns']))


class TestFilterUnique(basetest.BaseTestPycolator):
    command = 'filteruni'
    infilename = 'percolator_out.xml'
    suffix = '_filtuniq.xml'
    # FIXME other scores than svm
    # FIXME make sure BEST psm is retained, not worst.
    # FIXME illegal scores handling
    # FIXME PSM peptide reffing

    def test_filter_uniques(self):
        """Checks if resultpeps gets uniques, and also that input peptides
        were not unique to start with."""
        self.run_pycolator(['-s', 'svm'])
        result = self.read_percolator_out(self.resultfn)
        origin = self.read_percolator_out(self.infile)
        resultpeps = self.get_element_ids(result['peptides'],
                                          'peptide_id', result['ns'])
        originpeps = self.get_element_ids(origin['peptides'],
                                          'peptide_id', origin['ns'])
        self.assertEqual(len({x for x in resultpeps}), len(resultpeps))
        self.assertNotEqual(len({x for x in originpeps}), len(originpeps))


class TestFilterLength(basetest.BaseTestPycolator):
    command = 'filterlen'
    infilename = 'percolator_out.xml'
    suffix = '_filt_len.xml'
    # FIXME need to check maxlen minlen input?

    def strip_modifications(self, pep):
        return re.sub('\[UNIMOD:\d*\]', '', pep)

    def length_correct(self, seqs, minlen, maxlen):
        for seq in seqs:
            seq = self.strip_modifications(seq)
            self.assertGreaterEqual(len(seq), minlen)
            self.assertLessEqual(len(seq), maxlen)

    def all_peps_in_output(self, original_seqs, minlen, maxlen, testresult):
        for feat in original_seqs:
            seq = self.strip_modifications(feat)
            if len(seq) < minlen or len(seq) > maxlen:
                continue
            self.assertIn(feat, testresult)

    def test_filterlen(self):
        maxlen = 20
        minlen = 10
        self.run_pycolator(['--maxlen', str(maxlen), '--minlen', str(minlen)])
        result = self.get_psm_pep_ids_from_file(self.resultfn)
        origin = self.get_psm_pep_ids_from_file(self.infile)

        self.length_correct(result['peptide_ids'], minlen, maxlen)
        self.length_correct(result['psm_seqs'], minlen, maxlen)
        self.all_peps_in_output(origin['peptide_ids'], minlen, maxlen,
                                result['peptide_ids'])
        self.all_peps_in_output(origin['psm_seqs'], minlen, maxlen,
                                result['psm_seqs'])


class TestFilterNovel(basetest.LookupTestsPycolator):
    command = 'filterknown'
    infilename = 'percolator_out.xml'
    suffix = '_filtknown.xml'
    dbfn = 'known_peptide_lookup.sqlite'

    def test_noflags(self):
        self.dbpath = os.path.join(self.fixdir, self.dbfn)
        self.run_pycolator(['-b', self.dbpath])
        result = self.get_psm_pep_ids_from_file(self.resultfn)
        origin = self.get_psm_pep_ids_from_file(self.infile)
        self.assert_seqs_correct(origin['peptide_ids'], result['peptide_ids'])
        self.assert_seqs_correct(origin['psm_seqs'], result['psm_seqs'])

    def assert_seqs_correct(self, original_seqs, result_seqs):
        """Does the actual testing"""
        db = sqlite3.connect(self.dbpath)
        for oriseq in original_seqs:
            if self.seq_in_db(db, oriseq, '='):
                self.assertNotIn(oriseq, result_seqs)
            else:
                self.assertIn(oriseq, result_seqs)

    def test_ntermwildcards(self):
        self.dbpath = os.path.join(self.fixdir, self.dbfn)
        self.run_pycolator(['-b', self.dbpath, '--ntermwildcards'])
        assert 1 == 2


class TestTrypticLookup(basetest.LookupTestsPycolator):
    command = 'trypticlookup'
    infilename = 'proteins.fasta'
    suffix = '_lookup.sqlite'

    def query_db_assert(self, options=[], seqtype=None):
        with open(os.path.join(self.fixdir, 'peptides_trypsinized.yml')) as fp:
            tryp_sequences = yaml.load(fp)
        sequences = tryp_sequences['fully_tryptic']
        comparator = '='
        if seqtype == 'ntermfalloff':
            comparator = ' LIKE '
            sequences.extend(
                ['%{0}'.format(x) for x in tryp_sequences[seqtype]])
            sequences = [x[::-1] for x in sequences]
        elif seqtype is not None:
            sequences.extend(tryp_sequences[seqtype])
        self.run_pycolator(options)
        self.assertTrue(self.all_seqs_in_db(self.resultfn,
                                            sequences, comparator))

    def test_cutproline(self):
        self.query_db_assert(['--cutproline'], 'proline_cuts')

    def test_ntermwildcards(self):
        self.query_db_assert(['--ntermwildcards'], 'ntermfalloff')

    def test_noflags(self):
        self.query_db_assert()
