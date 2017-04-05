import os
from math import log

from tests.integration import basetests


class TestAddProtData(basetests.ProttableTest):
    command = 'proteindata'
    suffix = '_proteindata.txt'

    def setUp(self):
        super().setUp()
        self.dbfile = os.path.join(self.fixdir, 'mzidtsv_db.sqlite')

    def test_addprot(self):
        self.dbfile = os.path.join(self.fixdir, 'mzidtsv_db.sqlite')
        options = ['--setname', 'S1', '--dbfile', self.dbfile]
        self.run_command(options)
        self.check_protein_data('proteincentric')

    def test_genecentric(self):
        infn = 'prottable_genes.txt'
        self.infile = [os.path.join(self.fixdir, infn)]
        self.resultfn = os.path.join(self.workdir, infn + self.suffix)
        options = ['--setname', 'S1', '--dbfile', self.dbfile,
                   '--genecentric', 'genes']
        self.run_command(options)
        self.check_protein_data('genecentric')

    def test_assoccentric(self):
        infn = 'prottable_assoc.txt'
        self.infile = os.path.join(self.fixdir, infn)
        self.resultfn = os.path.join(self.workdir, infn + self.suffix)
        options = ['--setname', 'S1', '--dbfile', self.dbfile,
                   '--genecentric', 'assoc']
        self.run_command(options)
        self.check_protein_data('assoccentric')


class TestBuild(basetests.ProttableTest):
    command = 'build'
    suffix = ''
    infilename = 'built_protein_table.txt'

    def get_std_options(self):
        cmd = [self.executable, self.command, '-d', self.workdir]
        return cmd

    def test_proteincentric(self, cutoff=False):
        self.dbfile = os.path.join(self.fixdir, 'prottable_db.sqlite')
        options = ['--fdr', '--precursor', '--isobaric', '--probability',
                   '--pep', '--dbfile', self.dbfile]
        if cutoff:
            options.extend(['--mergecutoff', str(cutoff)])
        self.run_command(options)
        sql = ('SELECT p.protein_acc, bs.set_name, pf.fdr, pep.pep, pp.quant, '
               'ppr.probability FROM proteins AS p '
               'JOIN biosets AS bs '
               'JOIN protein_fdr AS pf USING(pacc_id) '
               'JOIN protein_pep AS pep USING(pacc_id) '
               'JOIN protein_precur_quanted AS pp USING(pacc_id) '
               'JOIN protein_probability AS ppr USING(pacc_id) '
               )

        self.check_build_values(sql, ['q-value', 'PEP', 'MS1 precursor area',
                                      'Protein error probability'],
                                'Protein accession', cutoff)
        sql = ('SELECT p.protein_acc, bs.set_name, pc.channel_name, '
               'pi.quantvalue, pi.amount_psms, pf.fdr FROM proteins AS p '
               'JOIN biosets AS bs '
               'JOIN protein_iso_quanted AS pi USING(pacc_id) '
               'JOIN protquant_channels AS pc USING(channel_id) '
               'JOIN protein_fdr AS pf WHERE pf.prottable_id=pc.prottable_id '
               'AND pf.pacc_id=p.pacc_id'
               )
        self.check_built_isobaric(sql, 'Protein accession', cutoff)
        self.check_protein_data('proteincentric')

    def test_mergecutoff(self):
        self.test_proteincentric(0.0001)

    def test_genecentric(self):
        self.dbfile = os.path.join(self.fixdir, 'prottable_gene_db.sqlite')
        options = ['--fdr', '--precursor', '--isobaric', '--probability',
                   '--pep', '--dbfile', self.dbfile, '--genecentric', 'genes']
        self.run_command(options)
        sql = ('SELECT g.gene_acc, bs.set_name, gf.fdr, pep.pep , gp.quant, '
               'gpr.probability '
               'FROM genes AS g '
               'JOIN biosets AS bs '
               'JOIN gene_fdr AS gf USING(gene_id) '
               'JOIN gene_pep AS pep USING(gene_id) '
               'JOIN gene_precur_quanted AS gp USING(gene_id) '
               'JOIN gene_probability AS gpr USING(gene_id) '
               )
        self.check_build_values(sql, ['q-value', 'PEP', 'MS1 precursor area',
                                      'Protein error probability'],
                                'Protein accession')
        sql = ('SELECT g.gene_acc, bs.set_name, pc.channel_name, '
               'pi.quantvalue, pi.amount_psms FROM genes AS g '
               'JOIN biosets AS bs '
               'JOIN gene_iso_quanted AS pi USING(gene_id) '
               'JOIN genequant_channels AS pc USING(channel_id) '
               )
        self.check_built_isobaric(sql, 'Protein accession')
        self.check_protein_data('genecentric')

    def test_assoccentric(self):
        self.dbfile = os.path.join(self.fixdir, 'prottable_assoc_db.sqlite')
        options = ['--fdr', '--precursor', '--isobaric', '--probability',
                   '--pep', '--dbfile', self.dbfile, '--genecentric', 'assoc']
        self.run_command(options)
        sql = ('SELECT ai.assoc_id, bs.set_name, gf.fdr, pep.pep , gp.quant, '
               'gpr.probability '
               'FROM associated_ids AS ai '
               'JOIN biosets AS bs '
               'JOIN assoc_fdr AS gf USING(gene_id) '
               'JOIN assoc_pep AS pep USING(gene_id) '
               'JOIN assoc_precur_quanted AS gp USING(gene_id) '
               'JOIN assoc_probability AS gpr USING(gene_id) '
               )
        self.check_build_values(sql, ['q-value', 'PEP', 'MS1 precursor area',
                                      'Protein error probability'],
                                'Protein accession')
        sql = ('SELECT ai.assoc_id, bs.set_name, pc.channel_name, '
               'pi.quantvalue, pi.amount_psms '
               'FROM associated_ids AS ai '
               'JOIN biosets AS bs '
               'JOIN assoc_iso_quanted AS pi USING(gene_id) '
               'JOIN genequant_channels AS pc USING(channel_id) '
               )
        self.check_built_isobaric(sql, 'Protein accession')
        self.check_protein_data('assoccentric')


class TestAddIsoquant(basetests.ProttableTest):
    command = 'addisoquant'
    suffix = '_added_isoq.txt'
    channels = ['tmt10plex_{}'.format(x) for x in ['126', '127N', '127C',
                                                   '128N', '128C', '129N',
                                                   '129C', '130N', '130C',
                                                   '131']]
    nopsms = ['{} - # quanted PSMs'.format(ch) for ch in channels]

    def test_isoquant(self):
        isotable = os.path.join(self.fixdir, 'prottable_isoquant.txt')
        options = ['--quantfile', isotable, '--isobquantcolpattern',
                   'tmt10plex', '--qaccpattern', 'accession']
        self.run_command(options)
        self.isoquant_check(isotable, 'Protein accession', self.channels,
                            self.nopsms)


class TestMS1Quant(basetests.ProttableTest):
    command = 'ms1quant'
    suffix = '_ms1q.txt'

    def setUp(self):
        super().setUp()
        self.psmfile = os.path.join(
            self.fixdir, 'mzidtsv_filtered_fr1-2_proteingrouped.txt')

    def check(self):
        top_ms1 = self.get_top_psms(self.psmfile, 'Peptide', 'MS1 area')
        top_ms1 = {prot: sum(sorted(ms1s.values(), reverse=True)[:3]) /
                   len(sorted(ms1s.values())[:3])
                   for prot, ms1s in top_ms1.items()}
        for protein in self.tsv_generator(self.resultfn):
            try:
                self.assertEqual(float(protein['MS1 precursor area']),
                                 top_ms1[protein['Protein accession']])
            except ValueError:
                self.assertNotIn(protein['Protein accession'], top_ms1)

    def test_ms1(self):
        options = ['--psmtable', self.psmfile]
        self.run_command(options)
        self.check()

    def test_with_protcol(self):
        options = ['--psmtable', self.psmfile, '--protcol', '14']
        self.run_command(options)
        self.check()


class TestProbability(basetests.ProttableTest):
    command = 'probability'
    suffix = '_protprob.txt'

    def check(self):
        top_psms = self.get_top_psms(self.pepfile, 'Peptide sequence', 'PEP',
                                     lowerbetter=True)
        for protein in self.tsv_generator(self.resultfn):
            protacc = protein['Protein accession']
            expected_prob = 1
            for prob in top_psms[protacc].values():
                expected_prob *= prob
            self.assertEqual(float(protein['Protein error probability']),
                             expected_prob)

    def setUp(self):
        super().setUp()
        self.pepfile = os.path.join(self.fixdir, 'peptable.txt')

    def test_proteincentric(self):
        options = ['--peptable', self.pepfile]
        self.run_command(options)
        self.check()

    def test_protcol(self):
        options = ['--peptable', self.pepfile, '--protcol', '15']
        self.run_command(options)
        self.check()


class TestEmpty(basetests.ProttableTest):
    command = 'emptytable'
    infilename = 'mzidtsv_filtered_fr1-2_proteingrouped.txt'
    suffix = '_prottable.txt'

    def check(self):
        expected_proteins = set()
        for psm in self.tsv_generator(self.infile[0]):
            prot = psm['Master protein(s)']
            if ';' in prot:
                continue
            expected_proteins.add(prot)
        res = [x['Protein accession'] for x in
               self.tsv_generator(self.resultfn)]
        self.assertFalse(expected_proteins.symmetric_difference(res))

    def test_proteingroup(self):
        self.run_command()
        self.check()

    def test_genecentric(self):
        options = ['--protcol', '14']
        self.run_command(options)
        self.check()


class TestBestpeptide(basetests.ProttableTest):
    command = 'bestpeptide'
    suffix = '_bestpep.tsv'

    def check(self):
        # FIXME what to do if a qvalue is at exactly 0? no log possible.
        # uglyfix here is to use PEP to avoid that scenario but the real
        # code has a nextbest qvalue for that.
        top_psms = self.get_top_psms(self.pepfile, 'Peptide sequence',
                                     'PEP', lowerbetter=True)
        qscores = {prot: -log(min([x for x in psms.values()]), 10)
                   for prot, psms in top_psms.items()}
        for protein in self.tsv_generator(self.resultfn):
            res = float(protein['Q-score best peptide'])
            self.assertEqual(qscores[protein['Protein accession']], res)

    def setUp(self):
        super().setUp()
        self.pepfile = os.path.join(self.fixdir, 'peptable.txt')

    def test_proteincentric(self):
        options = ['--logscore', '--scorecolpattern', '^PEP',
                   '--peptable', self.pepfile]
        self.run_command(options)
        self.check()

    def test_protcol(self):
        options = ['--logscore', '--scorecolpattern', '^PEP',
                   '--peptable', self.pepfile, '--protcol', '15']
        self.run_command(options)
        self.check()

    def test_protcol_pattern(self):
        options = ['--logscore', '--scorecolpattern', '^PEP',
                   '--peptable', self.pepfile, '--protcolpattern', 'Master']
        self.run_command(options)
        self.check()


class TestProtFDR(basetests.ProttableTest):
    command = 'protfdr'
    suffix = '_protfdr.txt'
    infilename = 'prottable_target.txt'

    def setUp(self):
        super().setUp()
        self.decoyfn = os.path.join(self.fixdir, 'prottable_decoy.txt')

    def test_qval(self):
        options = ['--decoyfn', self.decoyfn]
        self.run_command(options)
        expectedfn = os.path.join(self.fixdir, 'protfdr.txt')
        self.check_lines(expectedfn, self.resultfn)


class TestAssocPickFDR(basetests.ProttableTest):
    command = 'pickedfdr'
    suffix = '_pickedfdr.txt'
    infilename = 'prottable_assoc_target.txt'

    def test_result(self):
        self.infilename = 'prottable_assoc_target.txt'
        self.decoyfn = os.path.join(self.fixdir, 'prottable_assoc_decoy.txt')
        options = ['--decoyfn', self.decoyfn, '--picktype', 'result']
        self.run_command(options)
        expectedfn = os.path.join(self.fixdir, 'pickfdr_assoc.txt')
        self.check_lines(expectedfn, self.resultfn)


class TestFastaPickFDR(basetests.ProttableTest):
    command = 'pickedfdr'
    suffix = '_pickedfdr.txt'
    infilename = 'prottable_ensg_target.txt'

    def test_fasta(self):
        self.decoyfn = os.path.join(self.fixdir, 'prottable_ensg_decoy.txt')
        tfasta = os.path.join(self.fixdir, 'ensembl.fasta')
        dfasta = os.path.join(self.fixdir, 'decoy_ensembl.fasta')
        options = ['--decoyfn', self.decoyfn, '--targetfasta', tfasta,
                   '--decoyfasta', dfasta, '--picktype', 'fasta']
        self.run_command(options)
        expectedfn = os.path.join(self.fixdir, 'pickfdr_ensg.txt')
        self.check_lines(expectedfn, self.resultfn)


class TestAddFDR(basetests.ProttableTest):
    command = 'fdr'
    suffix = '_protfdr.txt'
    infilename = 'prottable_target.txt'

    def check_inbetween_fdr(self, score, fdrmap):
        def get_average(fdrmap, loscore, hiscore):
            low = [float(x) for x in fdrmap[str(loscore)]]
            hi = [float(x) for x in fdrmap[str(hiscore)]]
            return (str((low[0] + hi[0]) / 2), str((low[1] + hi[1]) / 2))

        for qv_score in sorted([float(x) for x in fdrmap.keys()]):
            if score > qv_score:
                lowerscore = str(qv_score)
            elif score < qv_score:
                return get_average(fdrmap, lowerscore, qv_score)

    def test(self):
        qvfn = os.path.join(self.fixdir, 'prottable_qvality.txt')
        options = ['-q', qvfn, '--scorecolpattern', '^Q-score']
        self.run_command(options)
        fdrmap = {float(q['Score']): (q['q-value'], q['PEP'])
                  for q in self.tsv_generator(qvfn)}
        for protein in self.tsv_generator(self.resultfn):
            score = round(float(protein['Q-score best peptide']), 5)
            try:
                self.assertEqual(protein['q-value'], fdrmap[score][0])
            except KeyError:
                fdr = self.check_inbetween_fdr(float(score), fdrmap)
                self.assertEqual(protein['q-value'], fdr[0])
                self.assertEqual(protein['PEP'], fdr[1])
            else:
                self.assertEqual(protein['PEP'], fdrmap[score][1])
