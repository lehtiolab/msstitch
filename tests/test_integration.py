import unittest
import subprocess
import os
import hashlib
from tempfile import mkdtemp
from lxml import etree


class TestSplitTD(unittest.TestCase):
    command = 'splittd'
    testdir = 'tests'
    fixdir = os.path.join(testdir, 'fixtures')
    outdir = os.path.join(testdir, 'test_output')
    infilename = 'percolator_out.xml'

    def setUp(self):
        self.infile = os.path.join(self.fixdir, self.infilename)
        os.makedirs(self.outdir, exist_ok=True)
        self.workdir = mkdtemp(dir=self.outdir)
        self.target = os.path.join(self.workdir,
                                   self.infilename + '_target.xml')
        self.decoy = os.path.join(self.workdir,
                                  self.infilename + '_decoy.xml')

    def tearDown(self):
        # remove self.workdir
        pass
    
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

    def run_pycolator(self, command, *options):
        cmd = ['./pycolator.py', '-c', command, '-i', self.infile,
               '-d', self.workdir]
        cmd.extend(options)
        subprocess.call(cmd)

    def test_split(self):
        """Tests that splitted files contain equal amount of PSMS
        when compared with expected output, and checks that each psm/peptide
        has correct 'decoy' attribute."""
        self.run_pycolator(self.command)
        target_expected = os.path.join(self.fixdir, 'splittd_target_out.xml')
        decoy_expected = os.path.join(self.fixdir, 'splittd_decoy_out.xml')
        target_exp_contents = self.read_percolator_out(target_expected)
        decoy_exp_contents = self.read_percolator_out(decoy_expected)
        target_contents = self.read_percolator_out(self.target)
        decoy_contents = self.read_percolator_out(self.decoy)
 
        self.assertEqual(len(target_contents['psms']), len(target_exp_contents['psms']))
        self.assertEqual(len(decoy_contents['peptides']), len(decoy_exp_contents['psms']))
        for feat in ['psms', 'peptides']:
            for el in target_contents[feat]:
                self.assertEqual(el.attrib['{%s}decoy' % target_contents['ns']], 'false')
            for el in decoy_contents[feat]:
                self.assertEqual(el.attrib['{%s}decoy' % decoy_contents['ns']], 'true')
