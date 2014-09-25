import basetests
from app.dataformats import mzidtsv as constants


class ProteinGroupTest(basetests.MzidTSVBaseTest):
    command = 'proteingroup'
    infilename = 'mzidtsv.txt'
    suffix = '_protgroups.txt'
    dbfn = 'proteingroup_lookup.sqlite'

    def run_and_analyze(self, options):
        self.run_command(options)
        result = self.get_proteingroups_per_peptide()
        self.do_asserting(result, self.expected)

    def get_proteingroups_per_peptide(self):
        with open(self.resultfn) as fp:
            header = next(fp).strip().split('\t')
            pep_ix = header.index(constants.HEADER_PEPTIDE)
            master_ix = header.index(constants.HEADER_MASTER_PROT)
            pgcontent_ix = header.index(constants.HEADER_PG_CONTENT)
            pgamount_ix = header.index(constants.HEADER_PG_AMOUNT_PROTEIN_HITS)
            for line in fp:
                line = line.strip().split('\t')
                print(line)
                yield {'peptide': line[pep_ix],
                       'master': line[master_ix],
                       'content': line[pgcontent_ix],
                       'amount': line[pgamount_ix],
                       }

    def do_asserting(self, result, expected, unrolled=False):
        pgkeys = ['peptide', 'master', 'content', 'amount']
        for line in expected:
            exp = {}
            for key, value in zip(pgkeys, line):
                exp[key] = value
            res = next(result)
            print(res)
            print(exp)
            self.assertEqual(res['peptide'], exp['peptide'])
            self.assertEqual(set(res['master'].split(';')),
                             set(exp['master'].split(';')))
            self.assertEqual(res['amount'], exp['amount'])
            rescontent = res['content'].split(';')
            expcontent = exp['content'].split(';')
            self.assertEqual(set(rescontent), set(expcontent))


class TestProteinGroupingConfidenceLevel(ProteinGroupTest):
    expected = [['ABCD', 'PB', 'PB,PA,PC', '2'],
                ['ABCD', 'PB', 'PB,PA,PC', '2'],
                ['HIJKL', 'PB;PD', 'PB,PA,PC;PD,PC,PE', '2;2'],
                ['MNOP', 'PD', 'PD,PC,PE', '1'],
                ['QRST', 'PD;PF', 'PD,PC,PE;PF,PG,PE', '2;2'],
                ['UVWX', 'PF', 'PF,PG,PE', '2'],
                ['YZAB', 'PF', 'PF,PG,PE', '2'],
                ['CDEF', 'PH', 'PH', '1'],
                ]

    def test_conflevel_higher_better(self):
        options = ['--confidence-col', '19', '--confidence-lvl', '0.01',
                   '--confidence-better', 'lower']
        self.run_and_analyze(options)

    def test_conflevel_lower_better(self):
        options = ['--confidence-col', '15', '--confidence-lvl', '3',
                   '--confidence-better', 'higher']
        self.run_and_analyze(options)


class TestProteinGroupingUnrolled(ProteinGroupTest):
    expected = [['K.ABCD.S', 'PB', 'PB,PA,PC', '2'],
                ['K.ABCD.S', 'PB', 'PB,PA,PC', '2'],
                ['ABCD', 'PB', 'PB,PA,PC', '2'],
                ['ABCD', 'PB', 'PB,PA,PC', '2'],
                ['HIJKL', 'PB;PD', 'PB,PA,PC;PD,PC,PE', '2;2'],
                ['HIJKL', 'PB;PD', 'PB,PA,PC;PD,PE,PC', '2;2'],
                ['HIJKL', 'PB;PD', 'PB,PA,PC;PD,PE,PC', '2;2'],
                ['.MNOP.', 'PD', 'PD,PC,PE', '1'],
                ['QRST', 'PD;PF', 'PD,PC,PE;PF,PG,PE', '2;2'],
                ['QRST', 'PD;PF', 'PD,PC,PE;PF,PG,PE', '2;2'],
                ['QRST', 'PD;PF', 'PD,PC,PE;PF,PG,PE', '2;2'],
                ['UVWX', 'PF', 'PF,PG,PE', '2'],
                ['UVWX', 'PF', 'PF,PG,PE', '2'],
                ['R.YZAB.L', 'PF', 'PF,PG,PE', '2'],
                ['R.YZAB.G', 'PF', 'PF,PG,PE', '2'],
                ['CDEF', 'PH', 'PH', '1'],
                ]

    infilename = 'mzidtsv_unrolled.txt'

    def test_unrolled(self):
        options = ['--confidence-col', '19', '--confidence-lvl', '0.01',
                   '--confidence-better', 'lower', '--unroll']
        self.run_and_analyze(options)
