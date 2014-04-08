import unittest
from app.lookups import sequences


class TestTrypsinize(unittest.TestCase):
    def runtest(self, seq, expected, proline=False):
        result = sequences.new_trypsinize(seq, proline)
        self.assertCountEqual(result, expected)

    def test_single_lysine(self):
        self.runtest('ABCAKABC', ['ABCAK', 'ABC'])

    def test_single_arginine(self):
        self.runtest('IAMAPEPKTIDE', ['IAMAPEPK', 'TIDE'])

    def test_proline_cuts(self):
        self.runtest('IAMARPEPKTIDE', ['IAMAR', 'PEPK', 'TIDE', 'IAMARPEPK'],
                     True)

    def test_proline_notcuts(self):
        self.runtest('IAMARPEPKTIDE', ['TIDE', 'IAMARPEPK'], False)

    def test_multitrypresidue(self):
        self.runtest('IAMAKKKEPT', ['IAMAK', 'K', 'K', 'EPT', 'IAMAKK',
                                    'IAMAKKK', 'KEPT', 'KKEPT', 'KK'])

    def test_multitrypresidue_with_proline(self):
        self.runtest('IAMAKKPKEPT', ['IAMAK', 'K', 'PK', 'KPK', 'EPT',
                                     'IAMAKK',
                                     'IAMAKKPK'], True)
