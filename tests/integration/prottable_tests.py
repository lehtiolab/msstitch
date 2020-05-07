import os
from math import log

from tests.integration import basetests


class TestBuild(basetests.ProttableTest):
    command = 'build'
    suffix = ''
    infilename = 'built_protein_table.txt'

    def get_std_options(self):
        cmd = [self.executable, self.command, '-d', self.workdir]
        return cmd

    def test_proteincentric(self, cutoff=False):
        self.dbfile = os.path.join(self.fixdir, 'proteins_db.sqlite')
        options = ['--fdr', '--precursor', '--isobaric', '--dbfile', self.dbfile]
        if cutoff:
            options.extend(['--mergecutoff', str(cutoff)])
        self.run_command(options)
        sql = ('SELECT p.protein_acc, bs.set_name, pf.fdr, pp.quant FROM proteins AS p '
               'JOIN biosets AS bs '
               'JOIN protein_fdr AS pf USING(pacc_id) '
               'JOIN protein_precur_quanted AS pp USING(pacc_id) '
               )
        self.check_build_values(sql, ['q-value', 'MS1 precursor area'],
                'Protein ID', cutoff)
        sql = ('SELECT p.protein_acc, bs.set_name, pc.channel_name, '
               'pi.quantvalue, pi.amount_psms, pf.fdr FROM proteins AS p '
               'JOIN biosets AS bs '
               'JOIN protein_iso_quanted AS pi USING(pacc_id) '
               'JOIN protquant_channels AS pc USING(channel_id) '
               'JOIN protein_fdr AS pf WHERE pf.prottable_id=pc.prottable_id '
               'AND pf.pacc_id=p.pacc_id'
               )
        self.check_built_isobaric(sql, 'Protein ID', cutoff=cutoff)
        self.check_protein_data('proteincentric')

    def test_mergecutoff(self):
        self.test_proteincentric(0.0001)

    def test_genecentric(self):
        self.dbfile = os.path.join(self.fixdir, 'genes_db.sqlite')
        options = ['--fdr', '--precursor', '--isobaric', 
                '--dbfile', self.dbfile, '--genecentric', 'genes']
        self.run_command(options)
        sql = ('SELECT g.gene_acc, bs.set_name, gf.fdr, gp.quant '
               'FROM genes AS g '
               'JOIN biosets AS bs '
               'JOIN gene_fdr AS gf USING(gene_id) '
               'JOIN gene_precur_quanted AS gp USING(gene_id) '
               )
        self.check_build_values(sql, ['q-value', 'MS1 precursor area'],
                'Gene ID')
        sql = ('SELECT g.gene_acc, bs.set_name, pc.channel_name, '
               'pi.quantvalue, pi.amount_psms FROM genes AS g '
               'JOIN biosets AS bs '
               'JOIN gene_iso_quanted AS pi USING(gene_id) '
               'JOIN genequant_channels AS pc USING(channel_id) '
               )
        self.check_built_isobaric(sql, 'Gene ID')
        self.check_protein_data('genecentric')

    def test_nopsmnrs(self):
        """Given a lookup with NO amount psm information, output a nice gene
        table without those columns"""
        self.dbfile = os.path.join(self.fixdir, 'proteins_nopsms_db.sqlite')
        options = ['--fdr', '--precursor', '--isobaric',
                '--dbfile', self.dbfile, '--genecentric', 'genes']
        self.run_command(options)
        sql = ('SELECT g.gene_acc, bs.set_name, gf.fdr, gp.quant '
               'FROM genes AS g '
               'JOIN biosets AS bs '
               'JOIN gene_fdr AS gf USING(gene_id) '
               'JOIN gene_precur_quanted AS gp USING(gene_id) '
               )
        self.check_build_values(sql, ['q-value', 'MS1 precursor area'],
                'Gene ID')
        sql = ('SELECT g.gene_acc, bs.set_name, pc.channel_name, '
               'pi.quantvalue, pi.amount_psms FROM genes AS g '
               'JOIN biosets AS bs '
               'JOIN gene_iso_quanted AS pi USING(gene_id) '
               'JOIN genequant_channels AS pc USING(channel_id) '
               )
        self.check_built_isobaric(sql, 'Gene ID', check_nr_psms=False)
        self.check_protein_data('genecentric')


    def test_assoccentric(self):
        self.dbfile = os.path.join(self.fixdir, 'assoc_db.sqlite')
        options = ['--fdr', '--precursor', '--isobaric',
                   '--dbfile', self.dbfile, '--genecentric', 'assoc']
        self.run_command(options)
        sql = ('SELECT ai.assoc_id, bs.set_name, gf.fdr, gp.quant '
               'FROM associated_ids AS ai '
               'JOIN biosets AS bs '
               'JOIN assoc_fdr AS gf USING(gene_id) '
               'JOIN assoc_precur_quanted AS gp USING(gene_id) '
               )
        self.check_build_values(sql, ['q-value', 'MS1 precursor area'], 
                'Gene Name')
        sql = ('SELECT ai.assoc_id, bs.set_name, pc.channel_name, '
               'pi.quantvalue, pi.amount_psms '
               'FROM associated_ids AS ai '
               'JOIN biosets AS bs '
               'JOIN assoc_iso_quanted AS pi USING(gene_id) '
               'JOIN genequant_channels AS pc USING(channel_id) '
               )
        self.check_built_isobaric(sql, 'Gene Name')
        self.check_protein_data('assoccentric')


class TestMS1Quant(basetests.ProttableTest):
    command = 'ms1quant'
    suffix = '_ms1q.txt'

    def setUp(self):
        super().setUp()
        self.psmfile = os.path.join(
            self.fixdir, 'target_pg.tsv')

    def check(self):
        top_ms1 = self.get_top_psms(self.psmfile, 'Peptide', 'MS1 area')
        top_ms1 = {prot: sum(sorted(ms1s.values(), reverse=True)[:3]) /
                   len(sorted(ms1s.values())[:3])
                   for prot, ms1s in top_ms1.items()}
        for protein in self.tsv_generator(self.resultfn):
            try:
                self.assertEqual(float(protein['MS1 precursor area']),
                                 top_ms1[protein['Protein ID']])
            except ValueError:
                self.assertNotIn(protein['Protein ID'], top_ms1)

    def test_ms1(self):
        options = ['--psmtable', self.psmfile]
        self.run_command(options)
        self.check()

    def test_with_protcol(self):
        options = ['--psmtable', self.psmfile, '--protcol', '14']
        self.run_command(options)
        self.check()


class TestBestpeptide(basetests.ProttableTest):
    command = 'bestpeptide'
    suffix = '_bestpep.tsv'

    def check(self, top_psms):
        # FIXME what to do if a qvalue is at exactly 0? no log possible.
        # uglyfix here is to use PEP to avoid that scenario but the real
        # code has a nextbest qvalue for that.
        qscores = {prot: -log(min([x for x in psms.values()]), 10)
                   for prot, psms in top_psms.items()}
        for protein in self.tsv_generator(self.resultfn):
            try:
                res = float(protein['Q-score best peptide'])
            except ValueError:
                self.assertEqual(protein['Q-score best peptide'], 'NA')
            else:
                self.assertEqual(qscores[protein['Protein ID']], res)

    def setUp(self):
        super().setUp()
        self.pepfile = os.path.join(self.fixdir, 'target_pep_quant.tsv')

    def test_proteincentric(self):
        options = ['--logscore', '--scorecolpattern', '^EValue',
                   '--peptable', self.pepfile]
        self.run_command(options)
        top_psms = self.get_top_psms(self.pepfile, 'Peptide sequence',
                                     'EValue', lowerbetter=True)
        self.check(top_psms)

    def test_only_NA(self):
        options = ['--logscore', '--scorecolpattern', '^FragMeth',
                   '--peptable', self.pepfile]
        self.run_command(options)
        top_psms = self.get_top_psms(self.pepfile, 'Peptide sequence',
                                     'FragMethod', lowerbetter=True)
        self.check(top_psms)

    def test_protcol(self):
        options = ['--logscore', '--scorecolpattern', '^EValue',
                   '--peptable', self.pepfile, '--protcol', '15']
        self.run_command(options)
        top_psms = self.get_top_psms(self.pepfile, 'Peptide sequence',
                                     'EValue', lowerbetter=True)
        self.check(top_psms)

    def test_protcol_pattern(self):
        options = ['--logscore', '--scorecolpattern', '^EValue',
                   '--peptable', self.pepfile, '--protcolpattern', 'Master']
        self.run_command(options)
        top_psms = self.get_top_psms(self.pepfile, 'Peptide sequence',
                                     'EValue', lowerbetter=True)
        self.check(top_psms)


class TestProtFDR(basetests.ProttableTest):
    command = 'protfdr'
    suffix = '_protfdr.txt'
    infilename = 'proteins_tbestpep'

    def setUp(self):
        super().setUp()
        self.decoyfn = os.path.join(self.fixdir, 'proteins_dbestpep')

    def test_qval(self):
        options = ['--decoyfn', self.decoyfn]
        self.run_command(options)
        expectedfn = os.path.join(self.fixdir, 'proteins_protfdr')
        self.check_lines(expectedfn, self.resultfn)


class TestAssocPickFDR(basetests.ProttableTest):
    command = 'pickedfdr'
    suffix = '_pickedfdr.txt'
    infilename = 'assoc_tbestpep'

    def test_result(self):
        self.decoyfn = os.path.join(self.fixdir, 'assoc_dbestpep')
        options = ['--decoyfn', self.decoyfn, '--picktype', 'result']
        self.run_command(options)
        expectedfn = os.path.join(self.fixdir, 'assoc_protfdr')
        self.check_lines(expectedfn, self.resultfn)


class TestFastaENSGPickFDR(basetests.ProttableTest):
    command = 'pickedfdr'
    suffix = '_pickedfdr.txt'
    infilename = 'genes_tbestpep'

    def test_ensg(self):
        self.decoyfn = os.path.join(self.fixdir, 'genes_dbestpep')
        tfasta = os.path.join(self.basefixdir, 'ens99_small.fasta')
        dfasta = os.path.join(self.fixdir, 'protrev_ens99_small.fasta')
        options = ['--decoyfn', self.decoyfn, '--targetfasta', tfasta,
                   '--decoyfasta', dfasta, '--picktype', 'ensg']
        self.run_command(options)
        expectedfn = os.path.join(self.fixdir, 'genes_protfdr')
        self.check_lines(expectedfn, self.resultfn)
    

class TestFastaGenenamePickFDR(basetests.ProttableTest):
    command = 'pickedfdr'
    suffix = '_pickedfdr.txt'
    infilename = 'assoc_tbestpep'

    def test_genename(self):
        self.decoyfn = os.path.join(self.fixdir, 'assoc_dbestpep')
        tfasta = os.path.join(self.basefixdir, 'ens99_small.fasta')
        dfasta = os.path.join(self.fixdir, 'protrev_ens99_small.fasta')
        options = ['--decoyfn', self.decoyfn, '--targetfasta', tfasta,
                   '--decoyfasta', dfasta, '--picktype', 'genename']
        self.run_command(options)
        expectedfn = os.path.join(self.fixdir, 'assoc_protfdr')
        self.check_lines(expectedfn, self.resultfn)
