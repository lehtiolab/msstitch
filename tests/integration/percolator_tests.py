import os
import sqlite3
from itertools import product

from tests.integration import basetests


class TestSplitProteinHeader(basetests.BaseTestPycolator):
    command = 'splitperco'
    suffix = ''

    def setUp(self):
        super().setUp()
        self.resultfn = None

    def test_splitprotein(self):
        protheaders = ['^ENSP000002', 'ENSP000002|^ENSP000003']
        self.suffix = ''
        options = ['--protheaders'] + protheaders
        self.run_command(options)
        ensp2_incl = os.path.join(self.workdir, self.infilename + '_h0.xml')
        ensp3_not_2 = os.path.join(self.workdir, self.infilename + '_h1.xml')
        first_expected = os.path.join(self.fixdir, 'perco.xml_h0.xml')
        second_expected = os.path.join(self.fixdir, 'perco.xml_h1.xml')
        incl_exp_contents = self.read_percolator_out(first_expected)
        excl_exp_contents = self.read_percolator_out(second_expected)
        incl_contents = self.read_percolator_out(ensp2_incl)
        excl_contents = self.read_percolator_out(ensp3_not_2)
        self.assertEqual(len(incl_contents['psms']),
                         len(incl_exp_contents['psms']))
        self.assertEqual(len(incl_contents['peptides']),
                         len(incl_exp_contents['peptides']))
        self.assertEqual(len(excl_contents['psms']),
                         len(excl_exp_contents['psms']))
        self.assertEqual(len(excl_contents['peptides']),
                         len(excl_exp_contents['peptides']))
        ns = self.get_xml_namespace(ensp2_incl)
        for feat in ['psms', 'peptides']:
            for el in incl_contents[feat]:
                self.assertIn(
                        'ENSP000002', [x.text[:10] for x in el.findall('{%s}protein_id' % ns['xmlns'])])
            for el in excl_contents[feat]:
                self.assertIn('ENSP000003', [x.text[:10] for x in el.findall('{%s}protein_id' % ns['xmlns'])])
                self.assertNotIn('ENSP000002', [x.text[:10] for x in el.findall('{%s}protein_id' % ns['xmlns'])])


class TestFilterKnown(basetests.BaseTestPycolator):
    command = 'filterperco'
    suffix = '_filtseq.xml'
    dbpath = 'seqs.db'

    def test_noflags(self):
        seqs = ['DSGTLQSQEAK']
        basetests.create_db(seqs)
        self.assert_seqs_correct(seqs)

    def test_ntermwildcards(self):
        seqs = ['DSGTLQSQEAK']
        falloff = 'LLLLLL'
        basetests.create_db([f'{falloff}{x}' for x in seqs], reverse=True)
        self.assert_seqs_correct(seqs, ['--insourcefrag', f'{len(falloff)+1}'])

    def test_deamidate_inmem(self):
        seqs = ['NSGTLQSQEAK']
        basetests.create_db(seqs)
        self.assert_seqs_correct(seqs, ['--deamidate', '--in-memory'], deamidate=True)

    def test_fullprotein(self):
        seqs = ['SEEAHRAEQLQDAEEEK']
        proteins = [f'XXXXXXXXX{x}XXXXXXXXX' for x in seqs]
        basetests.create_db(proteins, fullprotein=True, minlen=6)
        self.assert_seqs_correct(seqs, ['--fullprotein', '--minlen', '6'], fullprot=True)
        
    def deamidate(self, sequence):
        '''Deamidation is N->D, so any sequences found with D can have originally been N
        and should be checked with N in DB as well'''
        aa_possible = [(aa,) if aa != 'D' else ('D', 'N') for aa in sequence]
        return list(''.join(aa) for aa in product(*aa_possible))

    def assert_seqs_correct(self, seqs_filter, flags=[], deamidate=False, fullprot=False):
        """Does the actual testing"""
        options = ['--dbfile', self.dbpath]
        options.extend(flags)
        self.run_command(options)
        result = self.get_psm_pep_ids_from_file(self.resultfn)
        origin = self.get_psm_pep_ids_from_file(self.infile[0])
        seq_in_ori = {'peptide_ids': False, 'psm_seqs': False}
        for feattype in ['peptide_ids', 'psm_seqs']:
            original_seqs = origin[feattype]
            result_seqs = result[feattype]
            db = sqlite3.connect(self.dbpath)
            for oriseq in original_seqs:
                thisseq_filter = False
                seq_dbcheck = self.strip_modifications(oriseq)
                if deamidate:
                    testseqs = self.deamidate(seq_dbcheck)
                else:
                    testseqs = [seq_dbcheck]
                for seq in testseqs:
                    if seq in seqs_filter:
                        thisseq_filter = True
                        seq_in_ori[feattype] = True
                        self.assertNotIn(oriseq, result_seqs)
                if not thisseq_filter:
                    self.assertIn(oriseq, result_seqs)
        self.assertTrue(all(seq_in_ori.values()))



class TestDedupPeptides(basetests.BaseTestPycolator):
    command = 'dedupperco'
    suffix = '_dedup.xml'

    def test_do_no_duplicate_peps(self):
        options = []
        self.run_command(options)
        perco_input = self.get_psm_pep_ids_from_file(self.infile[0])
        perco_output = self.get_psm_pep_ids_from_file(self.resultfn)
        self.assertEqual(perco_input['psm_ids'], perco_output['psm_ids'])
        self.assertEqual(perco_input['peptide_ids'], perco_output['peptide_ids'])
        
    def test_remove_duplicate(self):
        self.infile = [os.path.join(self.basefixdir, 'perco_duplicate.xml')]
        options = []
        self.run_command(options)
        perco_input = self.get_psm_pep_ids_from_file(self.infile[0])
        perco_output = self.get_psm_pep_ids_from_file(self.resultfn)
        self.assertEqual(perco_input['psm_ids'], perco_output['psm_ids'])
        self.assertTrue(len(perco_input['peptide_ids']) > len(perco_output['peptide_ids']))
        self.assertEqual(set(perco_input['peptide_ids']), set(perco_output['peptide_ids']))

    def test_remove_duplicate_psms_too(self):
        self.infile = [os.path.join(self.basefixdir, 'perco_duplicate.xml')]
        options = ['--includepsms']
        self.run_command(options)
        perco_input = self.get_psm_pep_ids_from_file(self.infile[0])
        perco_output = self.get_psm_pep_ids_from_file(self.resultfn)
        self.assertEqual(len(perco_input['psm_ids']) - 1, len(perco_output['psm_ids']))
        self.assertEqual(len(perco_input['psm_seqs']) - 1, len(perco_output['psm_seqs']))
        self.assertEqual(set(perco_input['psm_ids']), set(perco_output['psm_ids']))
        self.assertEqual(set(perco_input['psm_seqs']), set(perco_output['psm_seqs']))
        self.assertEqual(len(perco_input['peptide_ids']) - 1, len(perco_output['peptide_ids']))
        self.assertEqual(set(perco_input['peptide_ids']), set(perco_output['peptide_ids']))
