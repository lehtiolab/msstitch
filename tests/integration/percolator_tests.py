import os
import sqlite3
from itertools import product

from tests.integration.basetests import BaseTestPycolator


class TestSplitProteinHeader(BaseTestPycolator):
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


class TestFilterKnown(BaseTestPycolator):
    command = 'filterperco'
    suffix = '_filtseq.xml'
    dbfn = 'known_peptide_lookup.sqlite'
    reversed_dbfn = 'rev_known_peptide_lookup.sqlite'

    def test_noflags(self):
        self.dbpath = os.path.join(self.basefixdir, self.dbfn)
        self.assert_seqs_correct()

    def test_ntermwildcards(self):
        max_falloff = 12
        self.dbpath = os.path.join(self.basefixdir, self.reversed_dbfn)
        self.assert_seqs_correct(['--insourcefrag', str(max_falloff)],
                                 'ntermfalloff', max_falloff)

    def test_deamidate(self):
        self.dbpath = os.path.join(self.basefixdir, self.dbfn)
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
