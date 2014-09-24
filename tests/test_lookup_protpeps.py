import unittest
from hashlib import md5
from unittest.mock import Mock, patch
from app.lookups import protein_peptide as lookup

DB_STORE_CHUNK = lookup.DB_STORE_CHUNK


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
        with patch('app.lookups.protein_peptide.sqlite.ProteinPeptideDB', self.mockdb), patch('app.lookups.protein_peptide.tsvreader.generate_tsv_psms', return_value=[]):
            lookup.create_protein_pep_lookup('testfn', [])
        self.mockdb.create_ppdb.assert_called_with()
        self.mockdb.store_peptides_proteins.assert_called_once_with({})
        self.mockdb.index.assert_called_once_with()


class TestCreateLookupLineParsing(TestCreateLookup):
    def pepprot_generator(self, line, *args):
        return self.specfn, line['scannr'], line['pepid'], line['seq'], line['proteins']

    def get_lines_and_expected(self, unroll):
        def get_line(scannr, pepseq, protein, pepid):
            line = {'pepid': pepid,
                    'scannr': scannr,
                    'seq': pepseq,
                   }
            if type(protein) is str:
                 line['proteins'] = [protein]
            else:
                 line['proteins'] = protein
            return line

        lines = []
        expected = {}
        scannr = 1234
        for pepseq, proteins in self.pepprots.items():
            pepid = md5('{0}{1}'.format(self.specfn, scannr)
                        .encode('utf-8')).hexdigest()
            if unroll:
                for protein in proteins:
                    lines.append(get_line(scannr, pepseq, protein, pepid))
            else:
                lines.append(get_line(scannr, pepseq, proteins, pepid))
            expected[pepid] = {'scan_nr': scannr, 'specfn': self.specfn,
                               'seq': pepseq, 'proteins': proteins}
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
                'app.lookups.protein_peptide.sqlite.ProteinPeptideDB',
                self.mockdb), patch('app.lookups.protein_peptide.tsvreader.generate_tsv_psms', return_value=lines):
            lookup.create_protein_pep_lookup('nofn', unroll)
        for expected_entry in expected:
            self.mockdb.store_peptides_proteins.assert_any_call(expected_entry)

    def test_multiprotein_lines(self):
        self.do_asserts(unroll=False)

    def test_unrolled_singleprotein_lines(self):
        self.do_asserts(unroll=True)

    def test_unrolled_singleprotein_small_chunk(self):
        self.do_asserts(unroll=True, chunk_size=1)
