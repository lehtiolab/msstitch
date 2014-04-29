import unittest
import subprocess
import os
import shutil
import hashlib
from tempfile import mkdtemp
from lxml import etree


class BaseTestPycolator(unittest.TestCase):
    testdir = 'tests'
    fixdir = os.path.join(testdir, 'fixtures')
    outdir = os.path.join(testdir, 'test_output')

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

    def run_pycolator(self, options=[]):
        cmd = ['./pycolator.py', '-c', self.command, '-i', self.infile,
               '-d', self.workdir]
        cmd.extend(options)
        subprocess.call(cmd)

    def setUp(self):
        self.infile = os.path.join(self.fixdir, self.infilename)
        os.makedirs(self.outdir, exist_ok=True)
        self.workdir = mkdtemp(dir=self.outdir)
        self.resultfn = os.path.join(self.workdir,
                                     self.infilename + self.suffix)

    def tearDown(self):
        shutil.rmtree(self.workdir)


class TestSplitTD(BaseTestPycolator):
    command = 'splittd'
    infilename = 'percolator_out.xml'

    def md5_check(self, fn):
        # DEPRECATE? XML too many formatting issues
        m = hashlib.md5()
        with open(fn) as fp:
            while True:
                data = fp.read(8192).encode('utf-8')
                if not data:
                    break
                m.update(data)
        return m.hexdigest()

    def test_split(self):
        """Tests that splitted files contain equal amount of PSMS
        when compared with expected output, and checks that each psm/peptide
        has correct 'decoy' attribute."""
        self.target = os.path.join(self.workdir,
                                   self.infilename + '_target.xml')
        self.decoy = os.path.join(self.workdir,
                                  self.infilename + '_decoy.xml')
        self.run_pycolator()
        target_expected = os.path.join(self.fixdir, 'splittd_target_out.xml')
        decoy_expected = os.path.join(self.fixdir, 'splittd_decoy_out.xml')
        target_exp_contents = self.read_percolator_out(target_expected)
        decoy_exp_contents = self.read_percolator_out(decoy_expected)
        target_contents = self.read_percolator_out(self.target)
        decoy_contents = self.read_percolator_out(self.decoy)

        self.assertEqual(len(target_contents['psms']), len(target_exp_contents['psms']))
        self.assertEqual(len(target_contents['peptides']), len(target_exp_contents['peptides']))
        self.assertEqual(len(decoy_contents['psms']), len(decoy_exp_contents['psms']))
        self.assertEqual(len(decoy_contents['peptides']), len(decoy_exp_contents['peptides']))
        for feat in ['psms', 'peptides']:
            for el in target_contents[feat]:
                self.assertEqual(el.attrib['{%s}decoy' % target_contents['ns']], 'false')
            for el in decoy_contents[feat]:
                self.assertEqual(el.attrib['{%s}decoy' % decoy_contents['ns']], 'true')


class TestMerge(BaseTestPycolator):
    command = 'merge'
    infilename = 'splittd_target_out.xml'
    suffix = '_merged.xml'

    def test_merge(self):
        self.multifiles = [os.path.join(self.fixdir, 'splittd_decoy_out.xml')]
        options = ['--multifiles']
        options.extend(self.multifiles)
        self.run_pycolator(options)
        expected = self.read_percolator_out(os.path.join(self.fixdir,
                                                         'percolator_out.xml'))
        result = self.read_percolator_out(self.resultfn)
        self.assertEqual(len(result['psms']), len(result['peptides']))
        self.assertCountEqual(self.get_element_ids(expected['psms'], 'psm_id',
                                                   expected['ns']),
                              self.get_element_ids(result['psms'], 'psm_id',
                                                   result['ns']))
        self.assertCountEqual(self.get_element_ids(expected['peptides'],
                                                   'peptide_id',
                                                   expected['ns']),
                              self.get_element_ids(result['peptides'],
                                                   'peptide_id', result['ns']))


class TestFilterUnique(BaseTestPycolator):
    command = 'filteruni'
    infilename = 'percolator_out.xml'
    suffix = '_filtuniq.xml'
    # FIXME other scores than svm
    # FIXME illegal scores handling
    # FIXME PSM peptide reffing

    def test_filter_uniques(self):
        """Checks if resultpeps gets uniques, and also that input peptides
        were not unique to start with."""
        self.run_pycolator(['-s', 'svm'])
        result = self.read_percolator_out(self.resultfn)
        origin = self.read_percolator_out(self.infile)
        resultpeps = self.get_element_ids(result['peptides'],
                                          'peptide_id', result['ns'])
        originpeps = self.get_element_ids(origin['peptides'],
                                          'peptide_id', origin['ns'])
        self.assertEqual(len({x for x in resultpeps}), len(resultpeps))
        self.assertNotEqual(len({x for x in originpeps}), len(originpeps))
