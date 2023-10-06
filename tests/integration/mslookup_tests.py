from tests.integration import basetests

import os
import sqlite3
import json
from Bio import SeqIO
from lxml import etree


class SearchspaceLookup(basetests.BaseTest):
    suffix = ''
    infilename = 'proteins.fasta'

    def setUp(self):
        super().setUp()
        self.infile = os.path.join(self.basefixdir, self.infilename)


class TestTrypsinize(SearchspaceLookup):
    suffix = '_tryp.fa'
    command = 'trypsinize'

    def run_case(self, minlen, cutproline, miscleav, ntermloss, ignorestopcodons):
        options = ['-o', self.resultfn]
        if minlen:
            options.extend(['--minlen', str(minlen)])
        if cutproline:
            options.extend(['--cutproline'])
            seqtype = 'proline_cuts'
        if ntermloss and miscleav and ignorestopcodons:
            options.extend(['--nterm-meth-loss', '--miscleav', str(miscleav), '--ignore-stop-codons'])
            seqtype = 'fully_tryptic_minlen_ntermloss_ignorestop'
        elif ntermloss and miscleav:
            options.extend(['--nterm-meth-loss', '--miscleav', str(miscleav)])
            seqtype = 'fully_tryptic_ntermloss'
        elif miscleav:
            options.extend(['--miscleav', str(miscleav)])
            seqtype = 'miscleav'
        self.run_command(options)
        with open(os.path.join(self.basefixdir, 'peptides_trypsinized.json')) as fp:
            tryp_sequences = json.load(fp)
        foundseqs = {}
        for rec in SeqIO.parse(self.resultfn, 'fasta'):
            self.assertEqual(tryp_sequences[seqtype][str(rec.seq)], rec.id)
            foundseqs[str(rec.seq)] = rec.id
            if minlen:
                self.assertGreaterEqual(len(str(rec.seq)), minlen)
        for pep, prot in tryp_sequences[seqtype].items():
            if not minlen or len(pep) >= minlen:
                self.assertEqual(prot, foundseqs[pep])

    def test_fullytryptic_minlen_ntermloss_stopcodons(self):
        self.run_case(8, False, 1, True, True)

    def test_prolinecut(self):
        self.run_case(False, True, False, False, False)

    def test_miss_cleavage(self):
        self.run_case(False, False, 1, False, False)

    def test_nterm_miss_cleavage(self):
        '''Specific code path for this, so test it'''
        self.run_case(False, False, 1, True, False)


class TestDecoyFa(SearchspaceLookup):
    command = 'makedecoy'
    infilename = 'twoproteins.fasta'

    def run_check(self, options):
        self.resultfn = os.path.join(self.workdir, 'decoy.fa')
        options.extend(['-o', self.resultfn])
        self.run_command(options)

    def run_without_predb(self, options):
        self.run_check(options)

    def run_with_existing_db(self, options):
        self.copy_db_to_workdir('decoycheck.sqlite', 'decoycheck.sqlite')
        options.extend(['--dbfile', 'decoycheck.sqlite'])
        self.run_check(options)

    def check_seqs(self, checkfile, dbcheck=False, keep_maxscrambled=False):
        targetfa = SeqIO.index(self.infile[0], 'fasta')
        checkfa = SeqIO.index(os.path.join(self.fixdir, checkfile), 'fasta')
        resfa = SeqIO.index(self.resultfn, 'fasta')
        ttryps = {}
        for tid, t in targetfa.items():
            ttryps[tid] = [y for x in t.seq.split('K') for y in x.split('R')]
        allttryps = (y for x in ttryps.values() for y in x)
        # Are all seqs from result in expected?
        for seqid, seq in resfa.items():
            check = checkfa[seqid]
            trypseqs = [y for x in seq.seq.split('K') for y in x.split('R')]
            trypchecks = [y for x in check.seq.split('K') for y in x.split('R')]
            try:
                self.assertEqual(seq.seq, check.seq)
            except AssertionError:
                if dbcheck:
                    tseqs = ttryps[seq.id.replace('decoy_', '')]
                    for resseq, tseq in zip(trypseqs, tseqs):
                        if resseq in allttryps:
                            self.assertEqual(resseq[-1], tseq[-1])
                            self.assertEqual(len(resseq), len(tseq))
                    self.assertEqual(set(seq.seq), set(check.seq))
                    self.assertEqual(len(seq), len(check))
                else:
                    raise
        # Are all seqs from check in result?
        for seqid, seq in checkfa.items():
            self.assertEqual(seq.seq, resfa[seqid].seq)

    def test_tryprev_predb_keeptargets(self):
        '''Tests making decoy where target hits are validated against DB, NOT removed if
        matching target despite shuffling'''
        self.run_with_existing_db(['--scramble', 'tryp_rev', '--maxshuffle', '10', '--keep-target'])
        self.check_seqs('decoy_tryprev_tcheck_keeptarget_twoproteins.fasta', dbcheck=True, keep_maxscrambled=True)

    def test_tryprev_yesdb(self):
        '''Tests making decoy where target hits are validated against DB, removed if
        matching target despite shuffling'''
        self.run_without_predb(['--scramble', 'tryp_rev'])
        self.check_seqs('decoy_tryprev_targetcheck_twoproteins.fasta', dbcheck=True)

    def test_tryprev_yesdb_minlen(self):
        '''Tests making decoy where target hits are validated against DB, removed if
        matching target despite shuffling and also removed if len < 5'''
        self.run_without_predb(['--scramble', 'tryp_rev', '--minlen', '5'])
        self.check_seqs('decoy_tryprev_minlen_twoproteins.fasta', dbcheck=True)

    def test_tryprev_ignore_db(self):
        '''Tests decoy where target hits are kept and no DB is checked'''
        self.run_without_predb(['--scramble', 'tryp_rev', '--ignore-target-hits'])
        self.check_seqs('decoy_tryprev_twoproteins.fasta', dbcheck=True)

    def test_protrev(self):
        '''Tests making decoy using protein-reverse, no target checking is done here'''
        self.run_without_predb(['--scramble', 'prot_rev'])
        self.check_seqs('decoy_twoproteins.fasta')

    def test_pretryp_predb(self):
        '''Tests making decoy where sequences are trypsinized before hand, and supplying
        an existing SQLite DB for matching target. Persistent (despite shuffling) target-matching 
        sequences will be removed'''
        options = ['--scramble', 'tryp_rev', '--notrypsin']
        self.infilename = 'trypsinized_twoproteins.fasta'
        self.infile = [os.path.join(self.fixdir, self.infilename)]
        self.run_with_existing_db(options)
        self.check_seqs('decoy_tryprev_pretryp_twoproteins.fasta', True)

    def test_pretryp_predb_keeptarget(self):
        '''Tests making decoy where sequences are trypsinized before hand, and supplying
        an existing SQLite DB for matching target. Persistent (despite shuffling) target-matching 
        sequences will not be removed'''
        options = ['--scramble', 'tryp_rev', '--notrypsin', '--keep-target']
        self.infilename = 'trypsinized_twoproteins.fasta'
        self.infile = [os.path.join(self.fixdir, self.infilename)]
        self.run_with_existing_db(options)
        self.check_seqs('decoy_tryprev_pretryp_keeptarget_twoproteins.fasta', True)


class TestTrypticLookup(SearchspaceLookup):
    command = 'storeseq'

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
        with open(os.path.join(self.basefixdir, 'peptides_trypsinized.json')) as fp:
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
        self.copy_db_to_workdir('spectra_lookup.sqlite')
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
    command = 'storeseq'

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
        with open(os.path.join(self.basefixdir, 'allpeptides_proteins.json')) as fp:
            sequences = [x.replace('L', 'I') for x in json.load(fp)]
        options.extend(['--fullprotein', '--minlen', '6'])
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
        self.copy_db_to_workdir('spectra_lookup.sqlite')
        self.query_db_assert(options)


class SpectraLookup(basetests.MSLookupTest):
    command = 'storespectra'

    def setUp(self):
        super().setUp()
        self.infile = os.path.join(self.basefixdir, self.infilename)

    def get_std_options(self):
        return [self.executable, self.command, '--spectra', self.infile]

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
        for scannr, spec in self.get_spectra_mzml([self.infile], bsets, ionmob):
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
                rtval = get_cvparam_value(scan, 'scan start time', ns)[0]
                if ionmob:
                    ion = get_cvparam_value(scan, 'inverse reduced ion mobility', ns)[0]
                    rt = float(rtval) / 60
                else:
                    ion = get_cvparam_value(scan, 'ion injection time', ns)[0]
                    rt = float(rtval)
                charge = get_cvparam_value(precursor, 'charge state', ns)
                charge = charge[0] if len(charge) else False
                mz = get_cvparam_value(precursor, 'selected ion m/z', ns)[0]
                exp_data = {'fn': os.path.basename(infile), 'bs': bset,
                            'charge': int(charge),
                            'sid': '{}_{}'.format(mzfnid + 1, scansid),
                            'mz': float(mz), 'rt': round(rt, 12), 'ion_something': round(float(ion), 12)}
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


class TestDDATIMSSpectraLookup(SpectraLookup):
    infilename = 'few_spec_timstof.mzML'

    def test_spectra_newdb(self):
        setnames = ['Set1']
        options = ['--setnames']
        options.extend(setnames)
        self.run_command(options)
        self.check_spectra(setnames, ionmob=True)


class TestSpecQuantLookup(basetests.MSLookupTest):
    command = 'storequant'
    infilename = ''
    isoinfilename = 'few_spectra.consXML'
    krinfilename = 'few_spectra.kr'
    dinoinfilename = 'few_spectra.dino'
    base_db_fn = 'spectra_lookup.sqlite'

    def setUp(self):
        super().setUp()
        self.isoinfile = os.path.join(self.fixdir, self.isoinfilename)
        self.krfile = os.path.join(self.fixdir, self.krinfilename)
        self.dinofile = os.path.join(self.fixdir, self.dinoinfilename)
        self.fakespfn = os.path.join(self.workdir, 'few_spectra.mzML')
        self.fakespfn2 = os.path.join(self.workdir, 'set2.mzML')

    def get_std_options(self):
        return [self.executable, self.command]

    def get_quantch_map(self):
        qmap_xml = {}
        for ac, map_el in etree.iterparse(self.isoinfile, tag='map'):
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
        sql = ('SELECT iq.intensity, ic.channel_name, pif.pif FROM isobaric_quant '
               'AS iq JOIN isobaric_channels AS ic USING(channel_id) '
               'LEFT OUTER JOIN precursor_ion_fraction AS pif USING(spectra_id)')
        qch_map = self.get_quantch_map()
        dbtmt = self.get_values_from_db(self.resultfn, sql)
        for ac, xml_quant in etree.iterparse(self.isoinfile,
                                             tag='consensusElement'):
            for element in xml_quant.findall('.//element'):
                qval, qchan, pif = next(dbtmt)
                self.assertEqual(qch_map[element.attrib['map']], qchan)
                self.assertEqual(float(element.attrib['it']), qval)
            xmlpif = xml_quant.xpath('./UserParam[@name="precursor_purity"]')
            if len(xmlpif):
                self.assertEqual(float(xmlpif[0].attrib['value']), pif)


    def test_isoquant(self):
        options = ['--isobaric', self.isoinfile, '--spectra', self.fakespfn]
        self.run_command(options)
        self.check_quantmap()
        self.check_quantification()

    def test_dinosaur(self):
        options = ['--dinosaur', self.dinofile, self.dinofile, '--rttol', '5', 
                '--mztol', '20', '--mztoltype', 'ppm', '--spectra', self.fakespfn, 
                self.fakespfn2]
        self.run_command(options)
        self.check_ms1_feats_stored(self.dinofile, 'dino', 'sum')

    def test_dinosaur_apex(self):
        options = ['--dinosaur', self.dinofile, self.dinofile, '--rttol', '5',
                '--mztol', '20', '--apex', '--mztoltype', 'ppm', '--spectra', self.fakespfn,
                self.fakespfn2]
        self.run_command(options)
        self.check_ms1_feats_stored(self.dinofile, 'dino', 'apex')

    def test_kronik(self):
        options = ['--kronik', self.krfile, self.krfile, '--rttol', '5',
                '--mztol', '20', '--mztoltype', 'ppm', '--spectra',
                self.fakespfn, self.fakespfn2]
        self.run_command(options)
        self.check_ms1_feats_stored(self.krfile, 'kr', 'sum')

    def test_kronik_apex(self):
        options = ['--kronik', self.krfile, self.krfile, '--rttol', '5',
                '--mztol', '20', '--apex', '--mztoltype', 'ppm', '--spectra',
                self.fakespfn, self.fakespfn2]
        self.run_command(options)
        self.check_ms1_feats_stored(self.krfile, 'kr', 'apex')

    def check_ms1_feats_stored(self, ms1file, feattype, intkey):
        PROTON_MASS = 1.0072
        sql = ('SELECT count(*) FROM ms1_quant LEFT OUTER JOIN ms1_fwhm USING(feature_id)')
        recs = self.get_values_from_db(self.resultfn, sql)
        sql = ('SELECT * FROM ms1_quant LEFT OUTER JOIN ms1_fwhm USING(feature_id)')
        feats = {1: {}, 2: {}}
        for recid, fnid, rt, mz, charge, val, fwhm in self.get_values_from_db(
                self.resultfn, sql):
            try:
                feats[fnid][rt][mz][charge] = (val, fwhm)
            except KeyError:
                try:
                    feats[fnid][rt][mz] = {charge: (val, fwhm)}
                except KeyError:
                    feats[fnid][rt] = {mz: {charge: (val, fwhm)}}
        fields = {
                'dino': {
                    'charge': 'charge',
                    'rt': 'rtApex',
                    'intensity': {'sum': 'intensitySum', 
                        'apex': 'intensityApex'}[intkey]
                    },
                'kr': {
                    'charge': 'Charge',
                    'rt': 'Best RTime',
                    'intensity': {'sum': 'Summed Intensity', 
                        'apex': 'Best Intensity'}[intkey]
                    },
                }[feattype]
        with open(ms1file) as fp:
            header = next(fp).strip().split('\t')
            for line in fp:
                line = {k: v for k, v in zip(header, line.strip('\n').split('\t'))}
                charge = float(line[fields['charge']])
                if feattype == 'kr':
                    line['fwhm'] = None
                    mz = (float(line['Monoisotopic Mass']) + charge * PROTON_MASS) / charge
                else:
                    mz = float(line['mz'])
                resval = feats[1][round(float(line[fields['rt']]), 12)][mz][charge]
                fwhm = float(line['fwhm']) if line['fwhm'] else None
                self.assertEqual((float(line[fields['intensity']]), fwhm), resval)
        # Now check we have a feature for all scans
        sql = """SELECT m.scan_sid, ma.feature_id FROM mzml AS m
        LEFT OUTER JOIN ms1_align AS ma USING(spectra_id)"""
        exemptscans = {'dino': [10229], 'kr': [10148, 10229]}
        exemptscans = [f'controllerType=0 controllerNumber=1 scan={x}' for x 
                in exemptscans[feattype]]
        for scan, featid in self.get_values_from_db(self.resultfn, sql):
            # some scans do not have a aligned MS1 value
            if scan not in exemptscans:
                self.assertFalse(featid == None)
            else:
                self.assertTrue(featid == None)
