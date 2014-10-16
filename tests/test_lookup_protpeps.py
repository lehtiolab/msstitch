import unittest
from unittest.mock import Mock, patch
from app.lookups import protein_peptide as lookup
import app.dataformats.mzidtsv as mzidtsvdata


DB_STORE_CHUNK = lookup.DB_STORE_CHUNK


class TestCreateLookup(unittest.TestCase):
    def setUp(self):
        self.mockdb = Mock
        self.mockdb.create_pgdb = Mock()
        self.mockdb.store_peptides_proteins = Mock()
        self.mockdb.index_protein_peptides = Mock()
        self.specfn = 'testfn'
        self.pepprots = {'IAMAPEPTIDE': ['ENS1234', 'ENS5678'],
                         'TRYIDENTIFYME': ['ENS9101112']}
        self.score = 999


class TestCreateLookupSqliteCalls(TestCreateLookup):
    def test_calls_sqlite(self):
        with patch('app.lookups.protein_peptide.ProteinGroupDB', self.mockdb), patch('app.lookups.protein_peptide.tsvreader.generate_tsv_psms', return_value=[]):
            lookup.create_protein_pep_lookup('testfn', [], 1, 1, True)
        self.mockdb.create_pgdb.assert_called_with()
        self.mockdb.store_peptides_proteins.assert_called_once_with({})
        self.mockdb.index_protein_peptides.assert_called_once_with()


class TestCreateLookupLineParsing(TestCreateLookup):
    def pepprot_generator(self, line, *args):
        return self.specfn, line[mzidtsvdata.HEADER_SCANNR], line['seq'], line['score'], line['proteins']

    def get_lines_and_expected(self, unroll):
        def get_line(scannr, pepseq, protein):
            line = {mzidtsvdata.HEADER_SCANNR: scannr,
                    mzidtsvdata.HEADER_SPECFILE: self.specfn,
                    'seq': pepseq,
                    'score': self.score,
                   }
            if type(protein) is str:
                 line['proteins'] = [protein]
            else:
                 line['proteins'] = protein
            return line

        lines = []
        expected = {}
        scannr = 1234
        rownr = 0
        for pepseq, proteins in self.pepprots.items():
            psm_id = '{0}_{1}'.format(self.specfn, scannr)
            expected[psm_id] = {'rows': [],
                                'seq': pepseq, 'proteins': proteins,
                                'score': self.score}
            if unroll:
                for protein in proteins:
                    lines.append(get_line(scannr, pepseq, protein))
                    expected[psm_id]['rows'].append(rownr)
                    rownr += 1
            else:
                lines.append(get_line(scannr, pepseq, proteins))
                expected[psm_id]['rows'].append(rownr)
                rownr += 1
            scannr += 1
        return lines, [expected]

    def do_asserts(self, unroll, chunk_size=DB_STORE_CHUNK):
        lines, expected = self.get_lines_and_expected(unroll)
        if chunk_size != DB_STORE_CHUNK:
            expected = [{k: v} for k, v in expected[0].items()]
        with patch(
            'app.lookups.protein_peptide.DB_STORE_CHUNK', chunk_size), patch(
                'app.lookups.protein_peptide.tsvreader.'
                'get_pepproteins', self.pepprot_generator), patch(
                'app.lookups.protein_peptide.ProteinGroupDB',
                self.mockdb), patch('app.lookups.protein_peptide.tsvreader.generate_tsv_psms', 
                return_value=lines), patch('app.lookups.protein_peptide.conffilt.passes_filter', return_value=True):
            lookup.create_protein_pep_lookup('nofn', ['header', 'header2'], 1, True, unroll)
        for expected_entry in expected:
            self.mockdb.store_peptides_proteins.assert_any_call(expected_entry)

    def test_multiprotein_lines(self):
        self.do_asserts(unroll=False)

    def test_unrolled_singleprotein_lines(self):
        self.do_asserts(unroll=True)

    def test_unrolled_singleprotein_small_chunk(self):
        self.do_asserts(unroll=True, chunk_size=1)
