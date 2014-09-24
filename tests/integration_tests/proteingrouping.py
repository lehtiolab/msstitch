import basetests
from app.dataformats import mzidtsv as constants


class ProteinGroupTest(basetests.BaseTest):
    command = 'proteingroup'
    infilename = 'mzidtsv.txt'
    suffix = '_proteingrouped.txt'
    dbfn = 'proteingroup_lookup.sqlite'

    def run_and_analyze(self, options):
        self.run_command(options)
        result = self.get_proteingroups_per_peptide()
        self.do_asserting(result, self.expected)

    def get_proteingroups_per_peptide(self):
        with open(self.resultfn) as fp:
            header = next(fp)
            pep_ix = header.index(constants.HEADER_PEPTIDE)
            master_ix = header.index(constants.HEADER_MASTER_PROT)
            pgcontent_ix = header.index(constants.HEADER_PG_CONTENT)
            pgamount_ix = header.index(constants.HEADER_PG_AMOUNT_PROTEIN_HITS)
            out = []
            for line in fp:
                out.append({'peptide': line[pep_ix],
                            'master': line[master_ix],
                            'content': line[pgcontent_ix],
                            'amount': line[pgamount_ix],
                            })
            return out

    def do_asserting(self, result, expected, unrolled=False):
        pgkeys = ['peptide', 'master', 'content', 'amount']
        for line in expected:
            exp = {}
            for key, value in zip(pgkeys, line[-1]):
                exp[key] = value
            if unrolled:
                exp['unroll'] = line[-1]
            else:
                exp['unroll'] = 1
            for pepline in range(exp['unroll']):
                res = next(result)
                for psm in result.items():
                    for key in pgkeys:
                        self.assertEqual(res[key], exp[key])


class TestProteinGroupingConfidenceLevel(ProteinGroupTest):
    expected = [['ABCD', 'PB', 'PA,PB,PC', '2', 2],
                ['ABCD', 'PB', 'PA,PB,PC', '2', 2],
                ['HIJKL', 'PB;PD', 'PA,PB,PC;PC,PD,PE', '2;2', 3],
                ['MNOP', 'PD', 'PC,PD,PE', '2', 1],
                ['QRST', 'PD;PF', 'PC,PD,PE;PE,PF,PG', '2,2', 3],
                ['UVWX', 'PF', 'PE,PF,PG', '2', 2],
                ['YZAB', 'PF', 'PE,PF,PG', '2', 2],
                ['CDEF', 'PH', 'PH', '1', 1],
                ]

    def test_conflevel_higher_better(self):
        options = ['--confidence-col', '19', '--confidence-level', '0.01',
                   '--confidence-better', 'lower']
        self.run_and_analyze(options)

    def test_conflevel_lower_better(self):
        options = ['--confidence-col', '15', '--confidence-level', '2',
                   '--confidence-better', 'higher']
        self.run_and_analyze(options)


class TestProteinGroupingUnrolled(ProteinGroupTest):
    infilename = 'mzidtsv_unrolled.txt'

    def test_unrolled(self):
        options = ['--confidence-col', '19', '--confidence-level', '0.01',
                   '--confidence-better', 'lower', '--unroll']
        self.run_and_analyze(options)
