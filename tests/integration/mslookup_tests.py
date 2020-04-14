from tests.integration import basetests

import os
import sqlite3
import json
from Bio import SeqIO
from lxml import etree


class SearchspaceLookup(basetests.BaseTest):
    suffix = ''
    infilename = 'proteins.fasta'
    executable = 'msslookup'


class TestTrypsinize(SearchspaceLookup):
    suffix = '_tryp.fa'
    command = 'trypsinize'
    infilename = 'proteins.fasta'

    def run_case(self, minlen, cutproline, miscleav):
        options = ['-o', self.resultfn]
        seqtype = 'fully_tryptic'
        if minlen:
            options.extend(['--minlen', str(minlen)])
        if cutproline:
            options.extend(['--cutproline'])
            seqtype = 'proline_cuts'
        if miscleav:
            options.extend(['--miscleav', str(miscleav)])
            seqtype = 'miscleav'
        cmd = self.run_command(options)
        with open(os.path.join(self.fixdir, 'peptides_trypsinized.yml')) as fp:
            tryp_sequences = json.load(fp)
        for rec in SeqIO.parse(self.resultfn, 'fasta'):
            self.assertEqual(tryp_sequences[seqtype][str(rec.seq)], rec.id)
            if minlen:
                self.assertGreaterEqual(len(str(rec.seq)), minlen)


    def test_fullytryptic(self):
        self.run_case(8, False, False)

    def test_prolinecut(self):
        self.run_case(False, True, False)

    def test_miss_cleavage(self):
        self.run_case(False, False, 1)


class TestDecoyFa(SearchspaceLookup):
    command = 'makedecoy'
    infilename = 'twoproteins.fasta'

    def run_check(self, options):
        self.resultfn = os.path.join(self.workdir, 'decoy.fa')
        options.extend(['-o', self.resultfn])
        self.run_command(options)

    def run_without_db(self, options):
        self.run_check(options)

    def run_with_existing_db(self, options):
        self.copy_db_to_workdir('decoycheck.sqlite', 'decoycheck.sqlite')
        options.extend(['--dbfile', 'decoycheck.sqlite'])
        self.run_check(options)

    def check_seqs(self, checkfile, targetscrambling=False):
        checkfa = SeqIO.index(os.path.join(self.fixdir, checkfile), 'fasta')
        resfa = SeqIO.index(self.resultfn, 'fasta')
        for seqid, seq in resfa.items():
            try:
                self.assertEqual(seq.seq, checkfa[seqid].seq)
            except AssertionError:
                if targetscrambling:
                    self.assertEqual(seq.seq[-1], checkfa[seqid].seq[-1])
                    self.assertEqual(set(seq.seq), set(checkfa[seqid].seq))
                else:
                    raise
        for seqid, seq in checkfa.items():
            try:
                self.assertEqual(seq.seq, resfa[seqid].seq)
            except AssertionError:
                if targetscrambling:
                    self.assertEqual(seq.seq[-1], resfa[seqid].seq[-1])
                    self.assertEqual(set(seq.seq), set(resfa[seqid].seq))
                else:
                    raise

    def test_tryprev_predb(self):
        self.run_with_existing_db(['--scramble', 'tryp_rev', '--maxshuffle', '10'])
        self.check_seqs('decoy_tryprev_targetcheck.fasta', targetscrambling=True)

    def test_tryprev_yesdb(self):
        self.run_without_db(['--scramble', 'tryp_rev'])
        self.check_seqs('decoy_tryprev_targetcheck.fasta', targetscrambling=True)

    def test_tryprev_yesdb_minlen(self):
        self.run_without_db(['--scramble', 'tryp_rev', '--minlen', '5'])
        self.check_seqs('decoy_tryprev_check_minlen.fasta', targetscrambling=True)

    def test_tryprev_ignore_db(self):
        self.run_without_db(['--scramble', 'tryp_rev', '--ignore-target-hits'])
        self.check_seqs('decoy_tryprev.fasta', targetscrambling=True)

    def test_protrev(self):
        self.run_without_db(['--scramble', 'prot_rev'])
        self.check_seqs('decoy_twoproteins.fasta')


class TestDecoyFaPretryp(SearchspaceLookup):
    command = 'makedecoy'
    infilename = 'twoproteins_tryp.fa'

    def test_tryprev_predb_trypsinized(self):
        self.infilename = 'twoproteins_tryp.fa'
        self.resultfn = os.path.join(self.workdir, 'decoy.fa')
        options = ['--scramble', 'tryp_rev', '--notrypsin', '-o', self.resultfn]
        self.run_command(options)
        self.check_seqs('decoy_twoproteins_tryp.fa', True)

    def check_seqs(self, checkfile, targetscrambling=False):
        checkfa = SeqIO.index(os.path.join(self.fixdir, checkfile), 'fasta')
        resfa = SeqIO.index(self.resultfn, 'fasta')
        for seqid, seq in resfa.items():
            try:
                self.assertEqual(seq.seq, checkfa[seqid].seq)
            except AssertionError:
                # peptide may have been shuffled when in db
                if targetscrambling:
                    self.assertEqual(seq.seq[-1], checkfa[seqid].seq[-1])
                    self.assertEqual(set(seq.seq), set(checkfa[seqid].seq))
                else:
                    raise


class TestTrypticLookup(SearchspaceLookup):
    command = 'seqspace'

    def all_seqs_in_db(self, dbfn, sequences, seqtype):
        db = sqlite3.connect(dbfn)
        seqs_in_db = set()
        for seq in sequences:
            seqs_in_db.add(self.seq_in_db(db, seq, seqtype))
        db.close()
        return seqs_in_db == set([True])

    def query_db_assert(self, options=None, seqtype=None):
        if options is None:
            options = []
        with open(os.path.join(self.fixdir, 'peptides_trypsinized.yml')) as fp:
            tryp_sequences = json.load(fp)
        if seqtype is not None:
            sequences = tryp_sequences[seqtype]
        else:
            sequences = tryp_sequences['fully_tryptic']
        self.run_command(options)
        self.assertTrue(self.all_seqs_in_db(self.resultfn,
                                            sequences, seqtype))

    def run_without_db(self, options=None, seqtype=None):
        self.resultfn = os.path.join(self.workdir,
                                     'mslookup_db.sqlite')
        self.query_db_assert(options, seqtype)

    def run_with_existing_db(self, options=None, seqtype=None):
        if options is None:
            options = []
        self.resultfn = os.path.join(self.workdir, 'seqspace.db')
        options.extend(['--dbfile', self.resultfn])
        self.copy_db_to_workdir('mzidtsv_db.sqlite')
        self.query_db_assert(options, seqtype)

    def test_cutproline_nodb(self):
        self.run_without_db(['--cutproline'], 'proline_cuts')

    def test_cutproline_yesdb(self):
        self.run_with_existing_db(['--cutproline'], 'proline_cuts')

    def test_ntermwildcards_nodb(self):
        self.run_without_db(['--insourcefrag'], 'ntermfalloff')

    def test_ntermwildcards_yes_db(self):
        self.run_with_existing_db(['--insourcefrag'], 'ntermfalloff')

    def test_noflags_no_db(self):
        self.run_without_db()

    def test_noflags_yes_db(self):
        self.run_with_existing_db()


class TestWholeProteinSeqLookup(SearchspaceLookup):
    command = 'protspace'

    def all_seqs_in_db(self, dbfn, sequences, seqtype):
        db = sqlite3.connect(dbfn)
        seqs_in_db = set()
        sql = ('SELECT EXISTS(SELECT seq FROM protein_peptides WHERE '
               'seq=? LIMIT 1)')
        for seq in sequences:
            seqs_in_db.add(db.execute(sql, (seq,)).fetchone()[0] == 1)
        db.close()
        return seqs_in_db == set([True])

    def query_db_assert(self, options):
        with open(os.path.join(self.fixdir, 'allpeptides_proteins.yml')) as fp:
            sequences = [x.replace('L', 'I') for x in json.load(fp)]
        options.extend(['--minlen', '6'])
        self.run_command(options)
        self.assertTrue(self.all_seqs_in_db(self.resultfn, sequences,
                                            seqtype=False))

    def test_without_db(self, seqtype=None):
        self.resultfn = os.path.join(self.workdir,
                                     'mslookup_db.sqlite')
        self.query_db_assert([])

    def test_with_existing_db(self, seqtype=None):
        self.resultfn = os.path.join(self.workdir, 'seqspace.db')
        options = ['--dbfile', self.resultfn]
        self.copy_db_to_workdir('mzidtsv_db.sqlite')
        self.query_db_assert(options)


class SpectraLookup(basetests.MSLookupTest):
    command = 'spectra'
    def check_spectra(self, bsets, ionmob=False):
        sql = ('SELECT mf.mzmlfilename, bs.set_name, s.scan_sid, s.charge, '
               's.mz, s.retention_time, s.spectra_id '
               '{} '
               'FROM mzml AS s '
               '{} '
               'JOIN mzmlfiles AS mf USING(mzmlfile_id) '
               'JOIN biosets AS bs USING(set_id)')
        if ionmob:
            ionsql = (', imob.ion_mobility', ' JOIN ionmob AS imob USING(spectra_id) ')
        else:
            ionsql = (', iit.ion_injection_time', ' JOIN ioninjtime AS iit USING(spectra_id) ')
        sql = sql.format(*ionsql)
        specrecs = {}
        for rec in self.get_values_from_db(self.resultfn, sql):
            specrecs[rec[2]] = {'fn': rec[0], 'bs': rec[1], 'charge': rec[3],
                                'mz': rec[4], 'rt': rec[5], 
                                'sid': rec[6],
                                'ion_something': rec[7],
                                }
        for scannr, spec in self.get_spectra_mzml(self.infile, bsets, ionmob):
            self.assertEqual(spec, specrecs[scannr])

    def get_spectra_mzml(self, infiles, bsets, ionmob):
        def get_cvparam_value(parent, name, ns):
            return [x.attrib['value'] for x in
                    parent.findall('{%s}cvParam' % ns['xmlns'])
                    if x.attrib['name'] == name]

        def multifind(elements, parent, ns):
            element = elements.pop(-1)
            if len(elements) > 0:
                parent = multifind(elements, parent, ns)
            return parent.find('{%s}%s' % (ns['xmlns'], element))

        for mzfnid, (bset, infile) in enumerate(zip(bsets, infiles)):
            ns = self.get_xml_namespace(infile)
            for ac, spectrum in etree.iterparse(
                    infile, tag='{%s}spectrum' % ns['xmlns']):
                if get_cvparam_value(spectrum, 'ms level', ns)[0] != "2":
                    continue
                scansid = spectrum.attrib['id']
                scan = multifind(['scanList', 'scan'], spectrum, ns)
                precursor = multifind(['precursorList', 'precursor',
                                       'selectedIonList', 'selectedIon'],
                                      spectrum, ns)
                rt = get_cvparam_value(scan, 'scan start time', ns)[0]
                if ionmob:
                    ion = get_cvparam_value(scan, 'inverse reduced ion mobility', ns)[0]
                else:
                    ion = get_cvparam_value(scan, 'ion injection time', ns)[0]
                print(etree.tostring(precursor))
                charge = get_cvparam_value(precursor, 'charge state', ns)
                charge = charge[0] if len(charge) else False
                mz = get_cvparam_value(precursor, 'selected ion m/z', ns)[0]
                exp_data = {'fn': os.path.basename(infile), 'bs': bset,
                            'charge': int(charge),
                            'sid': '{}_{}'.format(mzfnid + 1, scansid),
                            'mz': float(mz), 'rt': round(float(rt), 12), 'ion_something': round(float(ion), 12)}
                if not charge:
                    import sys; sys.stderr.write('{}\n'.format(charge))
                    sys.stderr.write('{}\n'.format(exp_data))
                yield scansid, exp_data


class TestDDAThermoSpectraLookup(SpectraLookup):
    infilename = 'few_spectra.mzML'

    def test_spectra_newdb(self):
        setnames = ['Set1']
        options = ['--setnames']
        options.extend(setnames)
        self.run_command(options)
        self.check_spectra(setnames)

    def test_spectra_newdb_filename(self):
        setnames = ['Set1']
        outfn = 'testresult.sql'
        self.resultfn = os.path.join(self.workdir, outfn)
        options = ['-o', outfn, '--setnames']
        options.extend(setnames)
        self.run_command(options)
        self.check_spectra(setnames)

    def test_spectra_olddb(self):
        self.copy_db_to_workdir('bioset_lookup.sqlite')
        setnames = ['Set1']
        options = ['--setnames']
        options.extend(setnames)
        self.run_command(options)
        self.check_spectra(setnames)


class TestDDATIMSSpectraLookup(SpectraLookup):
    infilename = 'few_spec_timstof.mzML'

    def test_spectra_newdb(self):
        setnames = ['Set1']
        options = ['--setnames']
        options.extend(setnames)
        self.run_command(options)
        self.check_spectra(setnames, ionmob=True)


class TestPSMLookup(basetests.MSLookupTest):
    command = 'psms'
    base_db_fn = 'spectra_lookup.sqlite'
    infilename = 'target.tsv'
    """DB and PSM table contain:
    - ENSEMBL proteins
    - a Uniprot swiss protein
    - A self-annotated protein
    - A non-annotated (only peptide) proteins
    """

    def check_db_fasta(self, fasta, exp_proteins=None, desc=True):
        if exp_proteins is None:
            exp_proteins = {}
            for rec in SeqIO.parse(fasta, 'fasta'):
                rd = rec.description
                gene = 'NA'
                if 'gene_symbol:' in rd:
                    six = rd.index('gene_symbol:') + 12
                    gix = rd.index('gene:') + 5 
                    gene = rd[gix: rd.index(' ', gix)]
                    desc = rd[rd.index('description:') + 12:]
                elif 'GN=' in rd:
                    six = rd.index('GN=') + 3
                    desc = [x for x in rd.split() if '=' not in x][1:]
                elif 'msstitch_fake_gene' in rd:
                    # special case fake fasta record for non-standard gene 
                    six, desc = False, 'NA'
                    symbol = rd.split()[-1]
                elif 'msstitch_fake_onlypeptide' in rd:
                    # special fake fasta record for unannotated peptide
                    six, symbol, desc = False, 'NA', 'NA'
                exp_proteins[rec.id] = {
                        'seq': rec.seq,
                        'gene': gene,
                        'desc': desc,
                        'symbol': rd[six: rd.index(' ', six)] if six else 'NA',
                        }
        self.check_db_base(exp_proteins)
        sql = ('SELECT ps.protein_acc, ps.sequence, g.gene_acc, aid.assoc_id, pd.description'
               ' FROM genes AS g JOIN associated_ids AS aid USING(protein_acc)'
               ' JOIN protein_seq AS ps USING(protein_acc)'
               ' JOIN prot_desc AS pd USING(protein_acc)')
        if not desc:
            sql = ('SELECT ps.protein_acc, ps.sequence '
                   'FROM protein_seq AS ps '
                   'JOIN prot_desc AS pd USING(protein_acc)')
        for prot, seq, gene, aid, desc in self.get_values_from_db(self.resultfn,
                                                             sql):
            self.assertEqual(exp_proteins[prot]['seq'], seq)
            if desc:
                self.assertEqual(exp_proteins[prot]['gene'], gene)
                self.assertEqual(exp_proteins[prot]['symbol'], aid)
                self.assertEqual(exp_proteins[prot]['desc'], desc)

    def check_db_base(self, expected_proteins=None):
        expected_psms = self.get_expected_psms()
        if expected_proteins is None:
            expected_proteins = {x for y in expected_psms.values()
                                 for x in y['proteins']}
        protsql = 'SELECT * FROM proteins'
        for protein in self.get_values_from_db(self.resultfn, protsql):
            self.assertIn(protein[1], expected_proteins)
        psmsql = ('SELECT ps.sequence, p.score, pr.rownr '
                  'FROM psmrows AS pr JOIN psms AS p USING(psm_id) '
                  'JOIN peptide_sequences AS ps USING(pep_id)')
        for psm in self.get_values_from_db(self.resultfn, psmsql):
            expected_psm = (expected_psms[psm[2]]['seq'],
                            expected_psms[psm[2]]['score'])
            self.assertEqual((psm[0], psm[1]), expected_psm)
        ppsql = ('SELECT pp.protein_acc, pr.rownr FROM psmrows AS pr '
                 'JOIN protein_psm AS pp USING(psm_id)')
        for protpsm in self.get_values_from_db(self.resultfn, ppsql):
            self.assertIn(protpsm[0], expected_psms[protpsm[1]]['proteins'])

    def get_expected_psms(self):
        header = self.get_tsvheader(self.infile[0])
        prot_ix = header.index('Protein')
        seq_ix = header.index('Peptide')
        score_ix = header.index('MSGFScore')
        psms = {}
        for row, line in enumerate(self.get_all_lines(self.infile[0])):
            line = line.strip('\n').split('\t')
            psms[row] = {'proteins': [x.split('(pre')[0] for x in
                                      line[prot_ix].split(';')],
                         'seq': line[seq_ix],
                         'score': line[score_ix],
                         }
        return psms

    def test_no_fasta(self):
        options = ['--spectracol', '1']
        self.run_command(options)
        self.check_db_base()

    def test_with_fasta(self):
        fastafn = os.path.join(self.fixdir, 'ens99_small.fasta')
        options = ['--spectracol', '1', '--fasta', fastafn, 
                   '--fastadelim', 'pipe', '--genefield', '2']
        self.run_command(options)
        self.check_db_fasta(fastafn)


class TestIsoquantLookup(basetests.MSLookupTest):
    command = 'isoquant'
    infilename = 'tmt_quant.consxml'
    base_db_fn = 'spectra_lookup.sqlite'

    def get_quantch_map(self):
        qmap_xml = {}
        for ac, map_el in etree.iterparse(self.infile[0], tag='map'):
            channeldata = map_el.attrib
            qmap_xml[channeldata['id']] = channeldata['label']
        return qmap_xml

    def check_quantmap(self):
        qmap_xml = self.get_quantch_map()
        qmap_xml = [qmap_xml[x]for x in sorted(qmap_xml.keys())]
        sql = 'SELECT channel_name FROM isobaric_channels'
        for xml_channel, db_rec in zip(qmap_xml,
                                       self.get_values_from_db(self.resultfn,
                                                               sql)):
            self.assertEqual(xml_channel, db_rec[0])

    def check_quantification(self):
        sql = ('SELECT iq.intensity, ic.channel_name FROM isobaric_quant '
               'AS iq JOIN isobaric_channels AS ic USING(channel_id)')
        qch_map = self.get_quantch_map()
        dbtmt = self.get_values_from_db(self.resultfn, sql)
        for ac, xml_quant in etree.iterparse(self.infile[0],
                                             tag='consensusElement'):
            for element in xml_quant.findall('.//element'):
                qval, qchan = next(dbtmt)
                self.assertEqual(qch_map[element.attrib['map']], qchan)
                self.assertEqual(float(element.attrib['it']), qval)

    def test_isoquant(self):
        fakespectra = os.path.join(self.workdir, 'task_0_dataset_17694.dat')
        options = ['--spectra', fakespectra]
        self.run_command(options)
        self.check_quantmap()
        self.check_quantification()


class TestMS1KronikLookup(basetests.MSLookupTest):
    command = 'ms1quant'
    infilename = 'kronik.txt'
    base_db_fn = 'spectra_lookup.sqlite'

    def test_kronik(self):
        self.fakespfn = os.path.join(self.workdir, 'task_0_dataset_17694.dat')
        options = ['--quanttype', 'kronik', '--rttol', '5', '--mztol', '20',
                   '--mztoltype', 'ppm', '--spectra', self.fakespfn]
        self.run_command(options)
        self.check_feats_stored()

    def check_feats_stored(self):
        PROTON_MASS = 1.0072
        sql = ('SELECT * FROM ms1_quant')
        feats = {1: {}}
        for recid, fnid, rt, mz, charge, val in self.get_values_from_db(
                self.resultfn, sql):
            try:
                feats[fnid][rt][mz][charge] = val
            except KeyError:
                try:
                    feats[fnid][rt][mz] = {charge: val}
                except KeyError:
                    feats[fnid][rt] = {mz: {charge: val}}
        with open(self.infile[0]) as fp:
            next(fp).strip().split('\t')
            for line in fp:
                line = line.strip('\n').split('\t')
                charge = float(line[3])
                mz = (float(line[4]) + charge * PROTON_MASS) / charge
                resval = feats[1][float(line[10])][mz][charge]
                self.assertEqual(float(line[6]), resval)


class TestMS1OpenMSLookup(basetests.MSLookupTest):
    command = 'ms1quant'
    infilename = 'featfinder.featXML'
    base_db_fn = 'spectra_lookup.sqlite'

    def test_featfind(self):
        self.fail('Not implemented tests for storing featfinder data yet')


class TestProteingroupLookup(basetests.MSLookupTest):
    command = 'proteingroup'
    infilename = 'mzidtsv_filtered_fr1-2.txt'
    base_db_fn = 'mzidtsv_db_nopgroups.sqlite'

    def test_proteingroup(self):
        self.run_command()
        self.check_pgmasters()
        self.check_coverage()
        self.check_pg_content()
        self.check_psm_protein_group()

    def check(self, sql, keyfun, valfun):
        exp_file = os.path.join(self.fixdir, 'mzidtsv_db.sqlite')
        result = {keyfun(x): valfun(x) for x in
                  self.get_values_from_db(self.resultfn, sql)}
        expected = {keyfun(x): valfun(x) for x in
                    self.get_values_from_db(exp_file, sql)}
        for key, value in result.items():
            self.assertIn(key, expected.keys())
            self.assertEqual(value, expected[key])

    def check_coverage(self):
        sql = 'SELECT * FROM protein_coverage'
        self.check(sql, lambda x: x[0], lambda x: x[1])

    def check_psm_protein_group(self):
        sql = ('SELECT ppg.psm_id, pgm.protein_acc FROM psm_protein_groups '
               'AS ppg JOIN protein_group_master AS pgm USING(master_id)')
        self.check(sql, lambda x: x[0], lambda x: x[1])

    def check_pg_content(self):
        sql = ('SELECT pgm.protein_acc, pgc.protein_acc, pgc.peptide_count, '
               'pgc.psm_count, pgc.protein_score '
               'FROM protein_group_content AS pgc '
               'JOIN protein_group_master AS pgm USING(master_id) '
               'ORDER BY pgm.protein_acc, pgc.protein_acc')
        self.check(sql, lambda x: x[0], lambda x: x[1:])

    def check_pgmasters(self):
        sql = ('SELECT * FROM protein_group_master')
        self.check(sql, lambda x: x[1], lambda x: 1)


class TestProtPepTableLookup(basetests.MSLookupTest):
    def run_checks(self, options, sql, checkfeats):
        self.run_command(options)
        result = self.get_values_from_db(self.resultfn, sql)
        with open(self.infile[0]) as fp:
            header = next(fp).strip().split('\t')
            for res, exp in zip(result, fp):
                exp = {field: val for field, val in
                       zip(header, exp.strip('\n').split('\t'))}
                for i, key in enumerate(checkfeats):
                    self.check_float(key, exp, res[i])

    def check_float(self, key, exp, result):
        expected = self.get_float_or_na(exp[key])
        self.assertEqual(result, expected)

    def check_isobaric(self, sql):
        def do_check(result, expected, header):
            expected = {k: v for k, v in
                        zip(header, expected.strip('\n').split('\t'))}
            [self.assertIn(key, expected) for key in result.keys()]
            [self.check_float(key, expected, result[key]) for key in result]
        protein_acc, resultquant = None, {}
        with open(self.infile[0]) as fp:
            header = next(fp).strip().split('\t')
            for rec in self.get_values_from_db(self.resultfn, sql):
                if rec[0] != protein_acc and protein_acc is not None:
                    do_check(resultquant, next(fp), header)
                    resultquant = {}
                    protein_acc = rec[0]
                elif protein_acc is None:
                    protein_acc = rec[0]
                resultquant[rec[1]] = rec[2]


class TestPeptidetableLookup(TestProtPepTableLookup):
    command = 'peptides'
    infilename = 'peptable.txt'
    base_db_fn = 'mzidtsv_db.sqlite'

    def run_peptable(self, opts=None):
        options = ['--peptidecol', '13', '--setnames', 'S1',
                   '--isobquantcolpattern', 'tmt10plex', '--fdrcolpattern',
                   '^q-value', '--pepcolpattern', '^PEP', '--psmnrcolpattern',
                   'quanted', '--ms1quantcolpattern', 'area']
        if opts is not None:
            options.extend(opts)
        sql = ('SELECT ps.sequence, ppq.quant, pf.fdr, pp.pep '
               'FROM peptide_precur_quanted AS ppq '
               'JOIN peptide_sequences AS ps USING(pep_id) '
               'JOIN peptide_fdr AS pf USING(pep_id) '
               'JOIN peptide_pep AS pp USING(pep_id) '
               )
        checkfeats = ['Peptide sequence', 'MS1 area (highest of all PSMs)',
                      'q-value', 'PEP']
        self.run_checks(options, sql, checkfeats)
        iso_sql = ('SELECT p.sequence, pc.channel_name, piq.quantvalue '
                   'FROM peptide_iso_quanted AS piq '
                   'JOIN peptide_sequences AS p USING(pep_id) '
                   'JOIN pepquant_channels AS pc USING(channel_id)')
        self.check_isobaric(iso_sql)

    def test_genecentric(self):
        self.run_peptable(['--genecentric'])

    def test_proteincentric(self):
        self.run_peptable()


class TestProteintableLookup(TestProtPepTableLookup):
    command = 'proteins'
    infilename = 'full_prottable.txt'
    base_db_fn = 'mzidtsv_db.sqlite'

    def test_proteincentric(self):
        options = ['--setnames', 'S1',
                   '--isobquantcolpattern', 'itraq4plex', '--fdrcolpattern',
                   'q-value', '--pepcolpattern', 'PEP', '--psmnrcolpattern',
                   'quanted', '--ms1quantcolpattern', 'area',
                   '--probcolpattern', 'probability']
        sql = ('SELECT p.protein_acc, ppq.quant, pf.fdr, pp.pep, '
               'ppr.probability '
               'FROM protein_precur_quanted AS ppq '
               'JOIN proteins AS p USING(pacc_id) '
               'JOIN protein_fdr AS pf USING(pacc_id) '
               'JOIN protein_pep AS pp USING(pacc_id) '
               'JOIN protein_probability AS ppr USING(pacc_id)')
        checkfeats = ['Protein accession', 'MS1 precursor area', 'q-value',
                      'PEP', 'Protein error probability']
        self.run_checks(options, sql, checkfeats)
        iso_sql = ('SELECT p.protein_acc, pc.channel_name, piq.quantvalue '
                   'FROM protein_iso_quanted AS piq '
                   'JOIN proteins AS p USING(pacc_id) '
                   'JOIN protquant_channels AS pc USING(channel_id)')
        self.check_isobaric(iso_sql)

    def test_genecentric_no_psmnrs(self):
        options = ['--setnames', 'S1', '--genecentric', 'genes',
                   '--isobquantcolpattern', 'itraq4plex', '--fdrcolpattern',
                   'q-value', '--pepcolpattern', 'PEP', 
                   '--ms1quantcolpattern', 'area', '--protcol', '2',
                   '--probcolpattern', 'probability']
        sql = ('SELECT p.gene_acc, ppq.quant, pf.fdr, pp.pep, '
               'ppr.probability '
               'FROM gene_precur_quanted AS ppq '
               'JOIN genes AS p USING(gene_id) '
               'JOIN gene_fdr AS pf USING(gene_id) '
               'JOIN gene_pep AS pp USING(gene_id) '
               'JOIN gene_probability AS ppr USING(gene_id)')
        checkfeats = ['Gene', 'MS1 precursor area', 'q-value',
                      'PEP', 'Protein error probability']
        self.run_checks(options, sql, checkfeats)
        iso_sql = ('SELECT p.gene_acc, pc.channel_name, piq.quantvalue '
                   'FROM gene_iso_quanted AS piq '
                   'JOIN genes AS p USING(gene_id) '
                   'JOIN genequant_channels AS pc USING(channel_id)')
        self.check_isobaric(iso_sql)

    def test_assoccentric(self):
        options = ['--setnames', 'S1', '--genecentric', 'assoc',
                   '--isobquantcolpattern', 'itraq4plex', '--fdrcolpattern',
                   'q-value', '--pepcolpattern', 'PEP', '--psmnrcolpattern',
                   'quanted', '--ms1quantcolpattern', 'area', '--protcol', '3',
                   '--probcolpattern', 'probability']
        sql = ('SELECT p.assoc_id, ppq.quant, pf.fdr, pp.pep, '
               'ppr.probability '
               'FROM assoc_precur_quanted AS ppq '
               'JOIN associated_ids AS p USING(gene_id) '
               'JOIN assoc_fdr AS pf USING(gene_id) '
               'JOIN assoc_pep AS pp USING(gene_id) '
               'JOIN assoc_probability AS ppr USING(gene_id)')
        checkfeats = ['Associated gene ID', 'MS1 precursor area', 'q-value',
                      'PEP', 'Protein error probability']
        self.run_checks(options, sql, checkfeats)
        iso_sql = ('SELECT p.assoc_id, pc.channel_name, piq.quantvalue '
                   'FROM assoc_iso_quanted AS piq '
                   'JOIN associated_ids AS p USING(gene_id) '
                   'JOIN genequant_channels AS pc USING(channel_id)')
        self.check_isobaric(iso_sql)
