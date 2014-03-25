import unittest
from unittest.mock import patch, Mock
from app.preparation import openmz


def mock_quantmaps(self):
    return ('116', '117', '118')


def mock_findquants(self, fn, scannr):
    return (('116', 1024), ('117', 4352), ('118', 4543))


def mock_tsv_header(fn):
    yield 'this\tis\ta\ttest\theader'.split('\t')


def mock_specscanlines(fn):
    yield ['this', 'is', 'a', 'line'], ('specfile', '1234')


class TestQuantDBLookups(unittest.TestCase):
    def setUp(self):
        self.mockdb = Mock
        self.mockdb.lookup_quant = mock_findquants
        self.mockdb.get_all_quantmaps = mock_quantmaps
        self.mockreader = Mock
        self.mockreader.get_mzidtsv_lines_scannr_specfn = mock_specscanlines

    def test_get_quantheader(self):
        fn = 'test'
        with patch('app.preparation.openmz.sqlite.QuantDB', self.mockdb):
            gqh = openmz.get_quant_header
            assert gqh(fn) == ['116', '117', '118']

    def test_lookup_quant(self):
        result = openmz.lookup_quant('fn', '1234', self.mockdb())
        assert result == {'116': 1024,
                          '117': 4352,
                          '118': 4543,
                          }

    def test_generate_mzidtsv_quanted(self):
        with patch('app.preparation.openmz.sqlite.QuantDB', self.mockdb), patch('app.preparation.openmz.readers', self.mockreader):
            psms = openmz.generate_psms_quanted('fn', 'tsvfn',
                                                list(mock_quantmaps(self)))
            result = next(psms)
        line = next(mock_specscanlines('test'))[0] + [
            x[1] for x in mock_findquants(self, 'test', '1234')]

        self.assertEqual(result, line)


class Test(unittest.TestCase):
    def test_convert_quantdict(self):
        qdata = {'116': 1024, '117': 2048}
        allin = openmz.convert_quantdata_to_line(qdata, ['116', '117'])
        self.assertEqual(allin, [1024, 2048])
        notfound = openmz.convert_quantdata_to_line(qdata,
                                                    ['116', '117', '118'])
        self.assertEqual(notfound, [1024, 2048, 'NA'])

    def create_tsv_header(self):
        qhead = list(mock_quantmaps(self))
        with patch('app.preparation.openmz.readers.get_tsv_header',
                   mock_tsv_header):
            result = openmz.create_tsv_header_quant('test', qhead)
        desired = next(mock_tsv_header()) + qhead
        self.assertEqual(result, desired)
