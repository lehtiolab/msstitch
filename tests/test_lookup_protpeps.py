import unittest
from unittest.mock import Mock, patch
from app.lookups import protein_peptides as lookup


class TestCreateLookup(unittest.TestCase):
    def setUp(self):
        self.mockdb = Mock

    def test_calls_sqlite(self):
        with patch(lookup.tsvreader.get_mzidtsv_lines_scannr_specfn,
                   return_value=[]), patch(lookup.sqlite.ProteinPeptideDB,
                                           self.mockdb):
            lookup.create_protein_pep_lookup('testfn')
        self.mockdb.create_ppdb.assert_called_with()
        self.mockdb.store_peptides_proteins.assert_called_once_with({})
        self.mockdb.index.assert_called_once_with()
