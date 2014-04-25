import unittest
import subprocess
import os
import hashlib
from tempfile import mkdtemp


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

    def md5_check(self, fn):
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
        self.run_pycolator(self.command)
        target_expected = 'tests_fixtures/splittd_target_out.xml'
        decoy_expected = 'tests_fixtures/splittd_decoy_out.xml'
        self.assertEqual(self.md5_check(self.target),
                         self.md5_check(target_expected))
        self.assertEqual(self.md5_check(self.decoy),
                         self.md5_check(decoy_expected))
