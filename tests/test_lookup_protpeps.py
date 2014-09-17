import unittest
from unittest.mock import Mock, patch
from app.lookups import protein_peptides as lookup


class TestCreateLookup(unittest.TestCase):
    def setUp(self):
        self.mockdb = Mock
        self.mockdb.create_ppdb = Mock()
        self.mockdb.store_peptides_proteins = Mock()
        self.mockdb.index = Mock()
        self.specfn = 'testfn'
        self.pepprots = {'IAMAPEPTIDE': ['ENS1234', 'ENS5678'],
                         'TRYIDENTIFYME': ['ENS9101112']}

    def test_calls_sqlite(self):
        with patch('app.lookups.protein_peptides.tsvreader.get_mzidtsv_lines_scannr_specfn',
                   return_value=[]), patch('app.lookups.protein_peptides.sqlite.ProteinPeptideDB',
                                           self.mockdb):
            lookup.create_protein_pep_lookup('testfn')
        self.mockdb.create_ppdb.assert_called_with()
        self.mockdb.store_peptides_proteins.assert_called_once_with({})
        self.mockdb.index.assert_called_once_with()

    def test_multiprotein_lines(self):
        lines = []
        scannr = 1234
        for pepseq, proteins in self.pepprots:
            line = ['f{0}'.format(x) for x in range(30)]
            line[2] = scannr
            line[9] = pepseq
            line[10] = ';'.join(['{0}(pre=R,post=S)'.format(x)
                                 for x in proteins])
            lines.append(line)
            scannr += 1

        with patch(
            'app.lookups.protein_peptides.tsvreader.'
            'get_mzidtsv_lines_scannr_specfn',
            return_value=self.scanfn_generator(lines)), patch(
                'app.lookups.protein_peptides.sqlite.ProteinPeptideDB',
                self.mockdb):
            lookup.create_protein_pep_lookup('nofn')

        self.mockdb.store_peptides_proteins.assert_called_with()

    def scanfn_generator(self, lines):
        for line in lines:
            scan_nr = line[2]
            yield line, (scan_nr, self.specfn)

    def test_unrolled_singleprotein_lines(self):
        assert 1 == 2
