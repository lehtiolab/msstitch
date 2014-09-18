import unittest
from hashlib import md5
from app.readers import tsv as reader


class TestGetPeptideProteins(unittest.TestCase):
    def setUp(self):
        self.line = ['f{0}'.format(x) for x in range(30)]
        self.specfn = 'testfn'
        self.scannr = '1234'
        self.pepseq = 'IAMAPEPTIDE'
        self.proteins = ['ENS12345', 'ENS6789']

    def do_asserts(self, line, unroll, proteins):
        exp_pepid = md5('{0}{1}'.format(self.specfn, self.scannr)
                        .encode('utf-8')).hexdigest()
        res_pepid, res_pepseq, res_proteins = reader.get_peptide_proteins(
            line, self.specfn, self.scannr, unroll)
        self.assertEqual(exp_pepid, res_pepid)
        self.assertEqual(proteins, res_proteins)
        self.assertEqual(self.pepseq, res_pepseq)

    def test_unroll(self):
        self.line[reader.TSV_PEPTIDE_COL] = 'K.{0}.T'.format(self.pepseq)
        self.line[reader.TSV_PROTEIN_COL] = self.proteins[0]
        self.do_asserts(self.line, True, [self.proteins[0]])

    def test_non_unroll(self):
        self.line[reader.TSV_PEPTIDE_COL] = self.pepseq
        self.line[reader.TSV_PROTEIN_COL] = '{0}(pre=K, post=T);{1}(pre=R, post=S)'.format(
            self.proteins[0], self.proteins[1])
        self.do_asserts(self.line, False, self.proteins)
