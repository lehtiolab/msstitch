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
    fixdir = os.path.join(os.getcwd(), testdir, 'fixtures')
    basefixdir = os.path.join(os.getcwd(), testdir, 'base_fixtures')
    outdir = os.path.join(os.getcwd(), testdir, 'test_output')
    executable = 'msstitch'

    def setUp(self):
        self.infile = os.path.join(self.fixdir, self.infilename)
        os.makedirs(self.outdir, exist_ok=True)
        self.workdir = mkdtemp(dir=self.outdir)
        os.chdir(self.workdir)
        self.resultfn = os.path.join(self.workdir, self.infilename)

    def tearDown(self):
        pass
        #shutil.rmtree(self.workdir)

    def get_std_options(self):
        cmd = [self.executable, self.command, '-i']
        if type(self.infile) != list:
            self.infile = [self.infile]
        cmd.extend(self.infile)
        if self.resultfn is not None:
            cmd.extend(['-o', self.resultfn])
        elif self.resultfn is None:
            cmd.extend(['-d', self.workdir])
        return cmd

    def run_command(self, options=[], return_error=False):
        cmd = self.get_std_options()
        cmd.extend(options)
        check = not return_error
        try:
            complete = subprocess.run(cmd, capture_output=True, text=True, check=check)
        except subprocess.CalledProcessError as e:
            print('Failed to run command: ', cmd)
            print('Output of command:')
            print(e.stdout)
            print(e.stderr)
            raise
        else:
            return complete

    def get_values_from_db(self, dbfile, sql):
        db = sqlite3.connect(dbfile)
        return db.execute(sql)

    def seq_in_db(self, dbconn, seq, seqtype, max_falloff=False):
        seq = seq.replace('L', 'I')
        if seqtype == 'ntermfalloff':
            seq = '{0}%'.format(seq[::-1])
            sql = 'SELECT seqs FROM known_searchspace WHERE seqs LIKE ?'
            for match in dbconn.execute(sql, (seq,)):
                if match[0][:-max_falloff] in seq:
                    return True
        else:
            sql = ('SELECT EXISTS(SELECT seqs FROM known_searchspace WHERE '
                   'seqs=? LIMIT 1)')
            return dbconn.execute(sql, (seq,)).fetchone()[0] == 1

    def get_tsvheader(self, fn):
        with open(fn) as fp:
            return next(fp).strip('\n').split('\t')

    def tsv_generator(self, fn):
        header = self.get_tsvheader(fn)
        tsv = self.get_all_lines(fn)
        for line in tsv:
            yield {field: val for field, val in
                   zip(header, line.strip('\n').split('\t'))}

    def get_all_lines(self, fn):
        with open(fn) as fp:
            next(fp)
            for line in fp:
                yield line

    def get_xml_namespace(self, fn):
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

    def get_root_el(self, fn):
        rootgen = etree.iterparse(fn, events=('start',))
        root = next(rootgen)[1]
        for child in root.getchildren():
            root.remove(child)
        return root

    def copy_db_to_workdir(self, dbfn, dst=False):
        if not dst: 
            shutil.copy(os.path.join(self.fixdir, dbfn), self.resultfn)
        else:
            shutil.copy(os.path.join(self.fixdir, dbfn), dst)

    def get_float_or_na(self, value):
        try:
            return float(value)
        except ValueError:
            return value

    def check_lines(self, expected, result):
        with open(expected) as fp, open(result) as resultfp:
            for expline, resline in zip(fp, resultfp):
                self.assertEqual(expline, resline)

    def isoquant_check(self, expected_isotable, acc_field, channels, nopsms):
        isoquant = {}
        accession = self.get_tsvheader(expected_isotable)[0]
        for line in self.tsv_generator(expected_isotable):
            acc = line.pop(accession)
            isoquant[acc] = line
        for result in self.tsv_generator(self.resultfn):
            for ch in channels + nopsms:
                try:
                    resval = float(result[ch])
                    expval = float(isoquant[result[acc_field]][ch])
                except ValueError:
                    # NA found
                    self.assertEqual(result[ch], isoquant[result[acc_field]][ch])
                else:
                    self.assertEqual(resval, expval)


class BaseTestPycolator(BaseTest):
    infilename = 'perco.xml'

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
        ns = self.get_xml_namespace(fn)['xmlns']
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

    def get_subelements(self, elements, subel, ns):
        return [element.find('{%s}%s' % (ns, subel)) for element in elements]

    def strip_modifications(self, pep):
        return re.sub('\[UNIMOD:\d*\]', '', pep)


class MzidTSVBaseTest(BaseTest):
    infilename = 'target.tsv'
    dbfn = 'target_psms.sqlite'

    def setUp(self):
        super().setUp()
        self.dbfile = os.path.join(self.fixdir, self.dbfn)

    def rowify(self, records):
        row, rownr = [], 0
        for record in records:
            if record[0] == rownr:
                row.append(record)
            else:
                yield row
                row = [record]
                rownr += 1

    def check_results_sql(self, checkfields, expected_values):
        for resultvals, exp_vals in zip(self.get_values(checkfields),
                                        expected_values):
            for resultval, expectval in zip(resultvals, exp_vals):
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

    def get_values(self, checkfields, outfile=False):
        if not outfile:
            outfile = self.resultfn
        with open(outfile) as fp:
            header = next(fp).strip('\n').split('\t')
            fieldindices = [header.index(field) for field in checkfields]
            row = 0
            for line in fp:
                line = line.strip('\n').split('\t')
                if len(checkfields) > 1:
                    yield [(row, field, line[ix]) for field, ix in
                           zip(checkfields, fieldindices)]
                else:
                    yield [(row, line[ix]) for ix in fieldindices]
                row += 1


class MSLookupTest(BaseTest):
    base_db_fn = None
    suffix = ''

    def setUp(self):
        super().setUp()
        self.resultfn = os.path.join(self.workdir, 'mslookup_db.sqlite')
        if self.base_db_fn is not None:
            self.copy_db_to_workdir(self.base_db_fn)

    def run_command(self, options=None):
        if options is None:
            options = []
        if self.base_db_fn is not None:
            options.extend(['--dbfile', self.resultfn])
        super().run_command(options)


class ProttableTest(BaseTest):

    def setUp(self):
        super().setUp()
        self.psmfile = os.path.join(
            self.fixdir, 'target_pg.tsv')
        self.decoyfn = os.path.join(self.fixdir, 'decoy_peptides.tsv')

    def check_ms1(self, featkey, featout):
        top_ms1 = self.get_top_peps(self.infile[0], featkey, 'Peptide sequence', 'MS1 area (highest of all PSMs)')
        top_ms1 = {prot: sum(sorted(ms1s.values(), reverse=True)[:3]) /
                   len(sorted(ms1s.values())[:3])
                   for prot, ms1s in top_ms1.items()}
        for protein in self.tsv_generator(self.resultfn):
            try:
                self.assertEqual(float(protein['MS1 precursor area']),
                                 top_ms1[protein[featout]])
            except ValueError:
                self.assertNotIn(protein['Protein ID'], top_ms1)

    def dotest_proteintable(self, scorecolpat, featkey, featout, summarize_method='denoms'):
        if summarize_method == 'denoms':
            summ = ['--denompatterns', '126']
        elif summarize_method == 'sweep':
            summ = ['--mediansweep']
        elif summarize_method == 'intensity':
            summ = ['--medianintensity']
        options = ['--scorecolpattern', scorecolpat,
                '--logscore', '--decoyfn', self.decoyfn, '--ms1quant',
                '--isobquantcolpattern', 'plex',
                *summ, '--psmtable', self.psmfile
                ] + self.specialoptions
        self.run_command(options)
        self.check_ms1(featkey, featout)

    def get_top_peps(self, fn, featkey, pepkey, valuekey, lowerbetter=False):
        top_vals = {}
        for pep in self.tsv_generator(fn):
            prot = pep[featkey]
            seq = pep[pepkey]
            try:
                value = float(pep[valuekey])
            except ValueError:
                continue
            if prot in top_vals and seq in top_vals[prot]:
                if lowerbetter and value > top_vals[prot][seq]:
                    continue
                elif not lowerbetter and value < top_vals[prot][seq]:
                    continue
            try:
                top_vals[prot][seq] = value
            except KeyError:
                top_vals[prot] = {seq: value}
        return top_vals


class MergeTest(BaseTest):
    infilename = ''
    command = 'merge'
    dbfn = 'target_psms.sqlite'

    def run_command(self, options, nopsms=False):
        self.infile = os.path.join(self.fixdir, self.infilename)
        self.resultfn = os.path.join(self.workdir, self.infilename)
        if not nopsms:
            self.options.extend(['--psmnrcolpattern', 'quanted'])
        super().run_command(options)

    def setUp(self):
        super().setUp()
        self.dbfile = os.path.join(self.workdir, self.dbfn)
        self.copy_db_to_workdir(self.dbfn, self.dbfile)
        self.options = ['--setnames', 'Set1', 
                '--fdrcolpattern', 'q-value', '--ms1quantcolpattern', 
                'area', '--isobquantcolpattern', 'plex', 
                '--dbfile', self.dbfile,
                ]

    def check_build_values(self, sql, fields, accession, cutoff=False):
        expected = {}
        for rec in self.get_values_from_db(self.dbfile, sql):
            # skip multiple entries of e.g protein group per peptide
            # or gene per protein
            try:
                expected[rec[0]][rec[1]] = rec[2:]
            except KeyError:
                expected[rec[0]] = {rec[1]: rec[2:]}
        for line in self.tsv_generator(self.resultfn):
            acc = line[accession]
            for setname, pepvals in expected[acc].items():
                for val, field in zip(pepvals, fields):
                    self.assertEqual(str(val),
                                     line['{}_{}'.format(setname, field)])
            expected.pop(acc)
        self.check_exp_empty(expected, cutoff)

    def check_exp_empty(self, expected, cutoff):
        if cutoff:
            self.assertFalse(expected == {})
        else:
            self.assertTrue(expected == {})

    def check_built_isobaric(self, sql, accession, nopsms=False, cutoff=False):
        expected = {}
        for rec in self.get_values_from_db(self.dbfile, sql):
            am_psm = rec[4]
            try:
                expected[rec[0]][rec[1]][rec[2]] = [rec[3], am_psm]
            except KeyError:
                try:
                    expected[rec[0]][rec[1]] = {rec[2]: [rec[3], am_psm]}
                except KeyError:
                    expected[rec[0]] = {rec[1]: {rec[2]: [rec[3], am_psm]}}
            if cutoff:
                expected[rec[0]][rec[1]][rec[2]] = [rec[3], am_psm, rec[5]]
        for line in self.tsv_generator(self.resultfn):
            for setname, fields in expected[line[accession]].items():
                for field, exp_val in fields.items():
                    setfield = '{}_{}'.format(setname, field)
                    if cutoff and exp_val[2] > cutoff:
                        exp_val = ['NA', 'NA']
                    self.assertAlmostEqual(float(line[setfield]), exp_val[0])
                    if not nopsms:
                        nr_psms = line['{} - # quanted PSMs'.format(setfield)]
                        self.assertEqual(nr_psms, str(exp_val[1]))
            expected.pop(line[accession])
        self.check_exp_empty(expected, cutoff)

    def check_protein_data(self, centrictype, sql, psm_sql):
        centric = {'proteincentric': 'pc', 'genecentric': 'gc',
                   'assoccentric': 'ac'}[centrictype]
        expected = {rec[0]: [x if x else 'NA' for x in rec[1:]] for rec in
                    self.get_values_from_db(self.dbfile, sql)}
        pdatalup = {
                'pc': {'acc': 'Protein ID', 'fields': [('Gene ID', 0), ('Gene Name', 1), ('Coverage', 3)]},
                'gc': {'acc': 'Gene ID', 'fields': [('Gene Name', 0), ('Protein ID(s)', 1)]},
                'ac': {'acc': 'Gene Name', 'fields': [('Gene ID', 0), ('Protein ID(s)', 1)]},
                }
        for row in self.tsv_generator(self.resultfn):
            acc = row[pdatalup[centric]['acc']]
            for (field, ix) in pdatalup[centric]['fields']:
                self.assertEqual(set(str(row[field]).split(';')), set(str(expected[acc][ix]).split(',')))
            self.assertEqual(row['Description'], expected[acc][2])
        expected, unipeps = {}, {}
        for rec in self.get_values_from_db(self.dbfile, psm_sql):
            pacc = rec[0]
            try:
                unipeps[rec[2]].add(pacc)
            except KeyError:
                unipeps[rec[2]] = set([pacc])
            try:
                expected[pacc]['psms'].add(rec[1])
            except KeyError:
                expected[pacc] = {'psms': set([rec[1]]), 'pep': set([rec[2]]),
                                  'unipep': 0}
            else:
                expected[pacc]['pep'].add(rec[2])
        for pep, prot in unipeps.items():
            if len(prot) == 1:
                expected[prot.pop()]['unipep'] += 1
        for protein in self.tsv_generator(self.resultfn):
            pacc = protein[pdatalup[centric]['acc']]
            poolname = 'Set1'
            self.assertEqual(protein['{}_# Unique peptides'.format(poolname)],
                             str(expected[pacc]['unipep']))
            self.assertEqual(protein['{}_# Peptides'.format(poolname)],
                             str(len(expected[pacc]['pep'])))
            self.assertEqual(protein['{}_# PSMs'.format(poolname)],
                             str(len(expected[pacc]['psms'])))

