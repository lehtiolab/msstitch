import os
import sqlite3
from itertools import product

from tests.integration.basetests import BaseTestPycolator


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


class TestSplit(BaseTestPycolator):
    def setUp(self):
        super().setUp()
        self.resultfn = None

    def do_check(self, first_exp_fn, sec_exp_fn, first_result, second_result):
        first_expected = os.path.join(self.fixdir, first_exp_fn)
        second_expected = os.path.join(self.fixdir, sec_exp_fn)
        first_exp_contents = self.read_percolator_out(first_expected)
        second_exp_contents = self.read_percolator_out(second_expected)
        first_contents = self.read_percolator_out(first_result)
        second_contents = self.read_percolator_out(second_result)

        self.assertEqual(len(first_contents['psms']),
                         len(first_exp_contents['psms']))
        self.assertEqual(len(first_contents['peptides']),
                         len(first_exp_contents['peptides']))
        self.assertEqual(len(second_contents['psms']),
                         len(second_exp_contents['psms']))
        self.assertEqual(len(second_contents['peptides']),
                         len(second_exp_contents['peptides']))
        return first_contents, second_contents


class TestSplitProteinHeader(TestSplit):
    command = 'splitprotein'
    suffix = ''

    def test_splitprotein(self):
        protheaders = ['ENSP', 'random;decoy']
        self.suffix = ''
        options = ['--protheaders'] + protheaders
        self.run_command(options)
        ensp_result = os.path.join(self.workdir, self.infilename + '_h0.xml')
        fake_result = os.path.join(self.workdir, self.infilename + '_h1.xml')
        ensp_contents, dr_contents = self.do_check('splitprotein_h0.xml',
                                                   'splitprotein_h1.xml',
                                                   ensp_result, fake_result)
        ns = self.get_xml_namespace(ensp_result)
        for feat in ['psms', 'peptides']:
            for el in ensp_contents[feat]:
                self.assertEqual(
                    el.find('{%s}protein_id' % ns['xmlns']).text[:4], 'ENSP')
            for el in dr_contents[feat]:
                self.assertIn(el.find('{%s}protein_id' % ns['xmlns']).text[:6],
                              ['decoy_', 'random'])


class TestSplitTD(TestSplit):
    command = 'splittd'
    suffix = ''

    def test_splittd(self):
        """TestSplitTD:: test_splittd
        Tests that splitted files contain equal amount of PSMS
        when compared with expected output, and checks that each psm/peptide
        has correct 'decoy' attribute."""
        target_result = os.path.join(self.workdir,
                                     self.infilename + '_target.xml')
        decoy_result = os.path.join(self.workdir,
                                    self.infilename + '_decoy.xml')
        self.run_command()
        t_contents, d_contents = self.do_check('splittd_target_out.xml',
                                               'splittd_decoy_out.xml',
                                               target_result, decoy_result)
        for feat in ['psms', 'peptides']:
            for el in t_contents[feat]:
                self.assertEqual(
                    el.attrib['{%s}decoy' % t_contents['ns']], 'false')
            for el in d_contents[feat]:
                self.assertEqual(
                    el.attrib['{%s}decoy' % d_contents['ns']], 'true')


class TestMerge(BaseTestPycolator):
    command = 'merge'
    infilename = 'splittd_target_out.xml'
    suffix = '_merged.xml'

    def test_merge_onefile(self):
        self.infile = [self.infile]
        self.run_command()
        expected = self.read_percolator_out(self.infile[0])
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
            if len(seq) <= minlen or len(seq) >= maxlen:
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
    command = 'filterseq'
    suffix = '_filtseq.xml'
    dbfn = 'known_peptide_lookup.sqlite'
    reversed_dbfn = 'rev_known_peptide_lookup.sqlite'

    def test_noflags(self):
        self.dbpath = os.path.join(self.fixdir, self.dbfn)
        self.assert_seqs_correct()

    def test_ntermwildcards(self):
        max_falloff = 12
        self.dbpath = os.path.join(self.fixdir, self.reversed_dbfn)
        self.assert_seqs_correct(['--insourcefrag', str(max_falloff)],
                                 'ntermfalloff', max_falloff)

    def test_deamidate(self):
        self.dbpath = os.path.join(self.fixdir, self.dbfn)
        self.assert_seqs_correct(['--deamidate'], 'deamidate')

    def deamidate(self, sequence):
        aa_possible = [(aa,) if aa != 'D' else ('D', 'N') for aa in sequence]
        return list(''.join(aa) for aa in product(*aa_possible))

    def assert_seqs_correct(self, flags=[], seqtype=None, max_falloff=False):
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
                if seqtype == 'deamidate':
                    testseqs = self.deamidate(seq_dbcheck)
                else:
                    testseqs = [seq_dbcheck]
                for seq in testseqs:
                    if self.seq_in_db(db, seq, seqtype, max_falloff):
                        self.assertNotIn(oriseq, result_seqs)
                        break
                    else:
                        self.assertIn(oriseq, result_seqs)


class TestQvality(BaseTestPycolator):
    command = 'qvality'
    suffix = '_qvalityout.txt'

    def test_qvality(self):
        self.infilename = 'splittd_target_out.xml'
        self.decoyfn = 'splittd_decoy_out.xml'
        self.decoy = os.path.join(self.fixdir, self.decoyfn)
        #options = ['--decoyfn', self.decoy, '--feattype', 'blaja', ]
        self.fail('pycolator qvality integration testing not implemented yet')
