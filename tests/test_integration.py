import unittest
import subprocess
import os
import hashlib
from tempfile import mkdtemp


class TestSplitTD(unittest.TestCase):
    command = 'splittd'

    def setUp(self):
        self.infile = 'tests/fixtures/percolator_out.xml'
        os.makedirs('tests/test_output', exist_ok=True)
        self.workdir = mkdtemp(dir='tests/test_output')
        self.target = os.path.basename(self.infile) + '_target.xml'
        self.decoy = os.path.basename(self.infile) + '_decoy.xml'

    def tearDown(self):
        # remove self.workdir
        pass

    def md5_check(self, fn):
        with open(fn):
            m = hashlib.md5(fn.read())
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
