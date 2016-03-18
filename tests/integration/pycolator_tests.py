import os
from tests.integration.basetests import BaseTestPycolator
import sqlite3


class TestReassign(BaseTestPycolator):
    command = 'reassign'
    suffix = '_reassigned.xml'

    def test_reassign_psm(self):
        self.run_qvality('psm')

    def test_reassign_peptide(self):
        self.run_qvality('peptide')

    def run_qvality(self, feattype):
        def get_mapvals_interpolated(qpmap, svm):
            testsvm = float(svm)
            for map_svm in sorted([float(x) for x in qpmap], reverse=True):
                if testsvm < map_svm:
                    upperval = str(map_svm)
                elif testsvm > map_svm:
                    lowerval = str(map_svm)
                    return tuple([str(sum([float(qpmap[upperval][x]),
                                           float(qpmap[lowerval][x])]) / 2)
                                  for x in [0, 1]])
        qvalityfn = os.path.join(self.fixdir, 'qvality.txt')
        self.run_command(['--qvality', qvalityfn, '--feattype', feattype])
        result = self.read_percolator_out(self.resultfn)
        origin = self.read_percolator_out(self.infile[0])
        qvmap = {}
        with open(qvalityfn) as fp:
            next(fp)  # skip header
            for line in fp:
                line = line.strip('\n').split('\t')
                qvmap[line[0]] = (line[1], line[2])
        feat = '{}s'.format(feattype)
        id_id = '{}_id'.format(feattype)
        resultids = self.get_element_ids(result[feat], id_id, result['ns'])
        resultsvms = self.get_svms(result[feat], result['ns'])
        resultqvals = self.get_qvals(result[feat], result['ns'])
        resultpeps = self.get_peps(result[feat], result['ns'])
        originids = self.get_element_ids(origin[feat], id_id, origin['ns'])
        originsvms = self.get_svms(origin[feat], origin['ns'])
        for oid, rid, osvm, rsvm, rq, rpep in zip(originids, resultids,
                                                  originsvms, resultsvms,
                                                  resultqvals, resultpeps):
            self.assertEqual(oid, rid)
            self.assertEqual(osvm, rsvm)
            try:
                mapvals = qvmap[str(round(float(osvm), 5))]
            except KeyError:
                mapvals = get_mapvals_interpolated(qvmap, osvm)
            self.assertEqual(rq, mapvals[1])
            self.assertEqual(rpep, mapvals[0])


class TestSplitTD(BaseTestPycolator):
    command = 'splittd'
    suffix = ''

    def test_split(self):
        """Tests that splitted files contain equal amount of PSMS
        when compared with expected output, and checks that each psm/peptide
        has correct 'decoy' attribute."""
        self.target = os.path.join(self.workdir,
                                   self.infilename + '_target.xml')
        self.decoy = os.path.join(self.workdir,
                                  self.infilename + '_decoy.xml')
        self.run_command()
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


class TestMerge(BaseTestPycolator):
    command = 'merge'
    infilename = 'splittd_target_out.xml'
    suffix = '_merged.xml'

    def test_merge(self):
        #self.multifiles = [os.path.join(self.fixdir, 'splittd_decoy_out.xml')]
        #self.multifiles =
        #options = ['--multifiles']
        #options.extend(self.multifiles)
        self.infile = [self.infile]
        self.infile.append(os.path.join(self.fixdir, 'splittd_decoy_out.xml'))
        self.run_command()
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


class TestFilterUnique(BaseTestPycolator):
    command = 'filteruni'
    suffix = '_filtuni.xml'
    # FIXME other scores than svm
    # FIXME make sure BEST psm is retained, not worst.
    # FIXME illegal scores handling
    # FIXME PSM peptide reffing

    def test_filter_uniques(self):
        """Checks if resultpeps gets uniques, and also that input peptides
        were not unique to start with."""
        self.run_command(['--score', 'svm'])
        result = self.read_percolator_out(self.resultfn)
        origin = self.read_percolator_out(self.infile[0])
        resultpeps = self.get_element_ids(result['peptides'],
                                          'peptide_id', result['ns'])
        originpeps = self.get_element_ids(origin['peptides'],
                                          'peptide_id', origin['ns'])
        self.assertEqual(len({x for x in resultpeps}), len(resultpeps))
        self.assertNotEqual(len({x for x in originpeps}), len(originpeps))


class TestFilterLength(BaseTestPycolator):
    command = 'filterlen'
    suffix = '_filtlen.xml'
    # FIXME need to check maxlen minlen input?

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
        self.run_command(['--maxlen', str(maxlen), '--minlen', str(minlen)])
        result = self.get_psm_pep_ids_from_file(self.resultfn)
        origin = self.get_psm_pep_ids_from_file(self.infile[0])

        self.length_correct(result['peptide_ids'], minlen, maxlen)
        self.length_correct(result['psm_seqs'], minlen, maxlen)
        self.all_peps_in_output(origin['peptide_ids'], minlen, maxlen,
                                result['peptide_ids'])
        self.all_peps_in_output(origin['psm_seqs'], minlen, maxlen,
                                result['psm_seqs'])


class TestFilterKnown(BaseTestPycolator):
    command = 'filterknown'
    suffix = '_filtknown.xml'
    dbfn = 'known_peptide_lookup.sqlite'
    reversed_dbfn = 'rev_known_peptide_lookup.sqlite'

    def test_noflags(self):
        self.dbpath = os.path.join(self.fixdir, self.dbfn)
        self.assert_seqs_correct()

    def test_ntermwildcards(self):
        self.dbpath = os.path.join(self.fixdir, self.reversed_dbfn)
        self.assert_seqs_correct(['--ntermwildcards'], 'ntermfalloff')

    def assert_seqs_correct(self, flags=[], seqtype=None):
        """Does the actual testing"""
        options = ['--dbfile', self.dbpath]
        options.extend(flags)
        self.run_command(options)
        result = self.get_psm_pep_ids_from_file(self.resultfn)
        origin = self.get_psm_pep_ids_from_file(self.infile[0])
        for feattype in ['peptide_ids', 'psm_seqs']:
            original_seqs = origin[feattype]
            result_seqs = result[feattype]
            db = sqlite3.connect(self.dbpath)
            for oriseq in original_seqs:
                seq_dbcheck = self.strip_modifications(oriseq)
                if self.seq_in_db(db, seq_dbcheck, seqtype):
                    self.assertNotIn(oriseq, result_seqs)
                else:
                    self.assertIn(oriseq, result_seqs)


class TestQvality(BaseTestPycolator):
    command = 'qvality'
    suffix = '_qvalityout.txt'

    def test_qvality(self):
        self.infilename = 'splittd_target_out.xml'
        self.decoyfn = 'splittd_decoy_out.xml'
        self.decoy = os.path.join(self.fixdir, self.decoyfn)

        options = ['--decoyfn', self.decoy, '--feattype', 'blaja', ]#'--qoptions']
        #self.run_command(options)
        self.fail('pycolator qvality integration testing not implemented yet')


