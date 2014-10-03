import unittest
from hashlib import md5
from app.readers import tsv as reader
from app.dataformats import mzidtsv as mzidtsvdata

class TestGetPeptideProteins(unittest.TestCase):
    def setUp(self):
        self.specfn = 'testfn'
        self.scannr = '1234'
        self.pepseq = 'IAMAPEPTIDE'
        self.proteins = ['ENS12345', 'ENS6789']
        self.line = {'f{0}'.format(x): 0 for x in range(30)}
        self.line[mzidtsvdata.HEADER_SPECFILE] = self.specfn
        self.line[mzidtsvdata.HEADER_SCANNR] = self.scannr
        self.line[mzidtsvdata.HEADER_MSGFSCORE] = 123

    def do_asserts(self, line, unroll, proteins):
        exp_pepid = md5('{0}{1}'.format(self.specfn, self.scannr)
                        .encode('utf-8')).hexdigest()
        res_specfn, res_scan, res_pepid, res_pepseq, res_proteins = \
            reader.get_pepproteins(line, unroll)
        self.assertEqual(exp_pepid, res_pepid)
        self.assertEqual(proteins, res_proteins)
        self.assertEqual(self.pepseq, res_pepseq)

    def test_unroll(self):
        self.line[mzidtsvdata.HEADER_PEPTIDE] = 'K.{0}.T'.format(self.pepseq)
        self.line[mzidtsvdata.HEADER_PROTEIN] = self.proteins[0]
        self.do_asserts(self.line, True, [self.proteins[0]])

    def test_non_unroll(self):
        self.line[mzidtsvdata.HEADER_PEPTIDE] = self.pepseq
        self.line[mzidtsvdata.HEADER_PROTEIN] = '{0}(pre=K, post=T);{1}(pre=R, post=S)'.format(
            self.proteins[0], self.proteins[1])
        self.do_asserts(self.line, False, self.proteins)
