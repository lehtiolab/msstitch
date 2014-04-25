import unittest
import subprocess
import os
from tempfile import mkdtemp


class TestSplitTD(unittest.TestCase):
    command = 'splittd'

    def setUp(self):
        self.infile = 'tests/fixtures/percolator_out.xml'
        self.workdir = mkdtemp(dir='tests/test_output')
        self.target = os.path.basename(self.infile) + '_target.xml'
        self.decoy = os.path.basename(self.infile) + '_decoy.xml'

    def tearDown(self):
        # remove self.workdir
        pass

    def run_pycolator(self, command, *options):
        cmd = ['pycolator.py', '-c', command, '-i', self.infile,
               '-d', self.workdir]
        cmd.extend(options)
        subprocess.call(cmd)

    def test_split(self):
        print self.infile
        print self.workdir
        print self.target
        print self.decoy
        self.run_pycolator(self.command)
