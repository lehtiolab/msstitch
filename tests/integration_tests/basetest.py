import unittest
import subprocess
import os
import shutil
from lxml import etree
from tempfile import mkdtemp


class BaseTestPycolator(unittest.TestCase):
    testdir = 'tests'
    fixdir = os.path.join(testdir, 'fixtures')
    outdir = os.path.join(testdir, 'test_output')

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
        return [psm.find('{%s}peptide_seq' % ns).attrib['seq'] for psm in psms]

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


