import unittest
import subprocess
import os
import shutil
import sqlite3
import re
from lxml import etree
from tempfile import mkdtemp


class BaseTest(unittest.TestCase):
    testdir = 'tests'
    fixdir = os.path.join(testdir, 'fixtures')
    outdir = os.path.join(testdir, 'test_output')

    def setUp(self):
        self.infile = os.path.join(self.fixdir, self.infilename)
        os.makedirs(self.outdir, exist_ok=True)
        self.workdir = mkdtemp(dir=self.outdir)
        self.resultfn = os.path.join(self.workdir,
                                     self.infilename + self.suffix)

    def tearDown(self):
        shutil.rmtree(self.workdir)

    def run_command(self, options=[]):
        cmd = ['python3', '{0}'.format(self.executable), self.command,
               '-d', self.workdir, '-i']
        if type(self.infile) != list:
            self.infile = [self.infile]
        cmd.extend(self.infile)
        cmd.extend(options)
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError:
            print('Failed to run executable {}'.format(self.executable))
            raise
        return cmd

    def run_command_expect_error(self, options=[]):
        try:
            cmd = self.run_command(options)
        except subprocess.CalledProcessError:
            pass
        else:
            self.fail('Command {} should throw an error'.format(cmd))


class BaseTestPycolator(BaseTest):
    executable = 'pycolator.py'
    infilename = 'percolator_out.xml'

    def get_psm_pep_ids_from_file(self, fn):
        contents = self.read_percolator_out(fn)
        return {'psm_ids': self.get_element_ids(contents['psms'],
                                                'psm_id', contents['ns']),
                'peptide_ids': self.get_element_ids(contents['peptides'],
                                                    'peptide_id',
                                                    contents['ns']),
                'psm_seqs': self.get_psm_seqs(contents['psms'],
                                              contents['ns'])
                }

    def read_percolator_out(self, fn):
        ns = self.get_namespace(fn)['xmlns']
        contents = {'ns': ns, 'psms': [], 'peptides': []}
        xml = etree.iterparse(fn)
        for ac, el in xml:
            if el.tag == '{%s}psm' % ns:
                contents['psms'].append(el)
            elif el.tag == '{%s}peptide' % ns:
                contents['peptides'].append(el)
        return contents

    def get_element_ids(self, elements, id_attrib, ns):
        return [x.attrib['{%s}%s' % (ns, id_attrib)] for x in elements]

    def get_psm_seqs(self, psms, ns):
        return [pepseq.attrib['seq'] for pepseq in
                self.get_subelements(psms, 'peptide_seq', ns)]

    def get_svms(self, elements, ns):
        return [svm.text for svm in
                self.get_subelements(elements, 'svm_score', ns)]

    def get_qvals(self, elements, ns):
        return [qval.text for qval in
                self.get_subelements(elements, 'q_value', ns)]

    def get_peps(self, elements, ns):
        return [pep.text for pep in
                self.get_subelements(elements, 'pep', ns)]

    def get_subelements(self, elements, subel, ns):
        return [element.find('{%s}%s' % (ns, subel)) for element in elements]

    def get_root_el(self, fn):
        rootgen = etree.iterparse(fn, events=('start',))
        root = next(rootgen)[1]
        for child in root.getchildren():
            root.remove(child)
        return root

    def get_namespace(self, fn):
        root = self.get_root_el(fn)
        ns = {}
        for prefix in root.nsmap:
            separator = ':'
            nsprefix = prefix
            if prefix is None:
                nsprefix = ''
                separator = ''
            ns['xmlns{0}{1}'.format(separator, nsprefix)] = root.nsmap[prefix]
        return ns

    def strip_modifications(self, pep):
        return re.sub('\[UNIMOD:\d*\]', '', pep)


class LookupTestsPycolator(BaseTestPycolator):
    def seq_in_db(self, dbconn, seq, seqtype):
        comparator = '='
        if seqtype == 'ntermfalloff':
            comparator = ' LIKE '
            seq = '{0}%'.format(seq[::-1])
        seq = seq.replace('L', 'I')
        sql = ('SELECT EXISTS(SELECT seqs FROM known_searchspace WHERE '
               'seqs{0}? LIMIT 1)'.format(comparator))
        return dbconn.execute(sql, (seq,)).fetchone()[0] == 1

    def all_seqs_in_db(self, dbfn, sequences, seqtype):
        db = sqlite3.connect(dbfn)
        seqs_in_db = set()
        for seq in sequences:
            seqs_in_db.add(self.seq_in_db(db, seq, seqtype))
        db.close()
        return seqs_in_db == set([True])


class MzidTSVBaseTest(BaseTest):
    executable = 'mzidtsv.py'
    infilename = 'mzidtsv.txt'

    def setUp(self):
        super().setUp()
        self.dbfile = os.path.join(self.fixdir, 'mzidtsv_db.sqlite')

    def get_values_from_db(self, dbfile, sql):
        db = sqlite3.connect(dbfile)
        return db.execute(sql)

    def get_all_lines(self, fn):
        with open(fn) as fp:
            next(fp)
            for line in fp:
                yield line

    def check_results(self, checkfields, expected_values):
        for resultvals in self.get_values(checkfields):
            for resultval, expectval in zip(resultvals, expected_values):
                self.assertEqual([str(x) if x is not None else 'NA'
                                  for x in expectval],
                                 [str(x) for x in resultval])

    def process_dbvalues_both(self, dbfile, sql, channel_fields, permfieldnrs,
                              permfieldnames):
        dbvals = self.get_values_from_db(dbfile, sql)
        outresults, rownr, permvals = [], 0, []
        for record in dbvals:
            if record[0] != rownr:
                for outresult in outresults:
                    yield outresult
                for pfname, pval in zip(permfieldnames, permvals):
                    yield (rownr, pfname, pval)
                if channel_fields != []:
                    outresults = [tuple([record[0]] +
                                        [record[x] for x in channel_fields])]
                permvals = [record[nr] for nr in permfieldnrs]
                rownr += 1
            else:
                permvals = [record[nr] for nr in permfieldnrs]
                if channel_fields != []:
                    result = [record[0]] + [record[x] for x in channel_fields]
                    outresults.append(tuple(result))

    def get_values(self, checkfields):
        with open(self.resultfn) as fp:
            header = next(fp).strip('\n').split('\t')
            fieldindices = [header.index(field) for field in checkfields]
            row = 0
            for line in fp:
                line = line.strip('\n').split('\t')
                if len(checkfields) > 1:
                    yield [(row, field, line[ix]) for field, ix in
                           zip(checkfields, fieldindices)]
                else:
                    yield [(row, line[ix]) for field, ix in
                           zip(checkfields, fieldindices)]
                row += 1