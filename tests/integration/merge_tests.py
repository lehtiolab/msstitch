import os
import re
from math import log

from tests.integration import basetests


class TestPeptideMerge(basetests.MergeTest):
    infilename = 'target_peptides.tsv'

    def test_proteincentric_in_memory(self):
        self.options.extend(['--pepcolpattern', 'peptide PEP', '--in-memory'])
        self.run_command(self.options)
        sql = ('SELECT ps.sequence, p.psm_id, prot.protein_acc, pd.description, '
               'g.gene_acc, aid.assoc_id, pc.coverage, pseq.sequence '
               'FROM peptide_sequences AS ps '
               'JOIN psms AS p USING(pep_id) '
               'JOIN psm_protein_groups USING(psm_id) '
               'JOIN protein_group_master AS pm USING(master_id) '
               'JOIN proteins AS prot USING(pacc_id) '
               'JOIN protein_seq AS pseq USING(protein_acc) '
               'LEFT OUTER JOIN prot_desc AS pd USING(pacc_id) '
               'LEFT OUTER JOIN ensg_proteins AS egp USING(pacc_id) '
               'LEFT OUTER JOIN genes AS g USING(gene_id) '
               'LEFT OUTER JOIN genename_proteins AS gnp USING(pacc_id) '
               'LEFT OUTER JOIN associated_ids AS aid USING(gn_id) '
               'JOIN protein_coverage AS pc USING(protein_acc) ')
        self.check_iso_and_peptide_relations(sql, proteincentric=True)

    def test_genecentric_in_memory(self):
        self.options.append('--genecentric')
        self.options.extend(['--pepcolpattern', 'peptide PEP', '--in-memory'])
        self.options.extend(['--flrcolpattern', 'q-value'])
        self.run_command(self.options)
        sql = ('SELECT ps.sequence, p.psm_id, "NA", pd.description, '
               'g.gene_acc, aid.assoc_id, "NA", "NA" '
               'FROM peptide_sequences AS ps '
               'JOIN psms AS p USING(pep_id) '
               'JOIN protein_psm USING(psm_id) '
               'JOIN proteins AS prot USING(protein_acc) '
               'LEFT OUTER JOIN prot_desc AS pd USING(pacc_id) '
               'LEFT OUTER JOIN ensg_proteins AS egp USING(pacc_id) '
               'LEFT OUTER JOIN genes AS g USING(gene_id) '
               'LEFT OUTER JOIN genename_proteins AS gnp USING(pacc_id) '
               'LEFT OUTER JOIN associated_ids AS aid USING(gn_id) '
               'JOIN protein_coverage AS pc USING(protein_acc) '
               )
        self.check_iso_and_peptide_relations(sql, flr=True)

    def test_nogroups(self):
        self.options.extend(['--pepcolpattern', 'peptide PEP'])
        self.options.append('--no-group-annotation')
        self.run_command(self.options)
        sql = ('SELECT ps.sequence, p.psm_id, prot.protein_acc, "NA", "NA", "NA", "NA", pseq.sequence '
               'FROM peptide_sequences AS ps '
               'JOIN psms AS p USING(pep_id) '
               'JOIN protein_psm USING(psm_id) '
               'JOIN proteins AS prot USING(protein_acc) '
               'JOIN protein_seq AS pseq USING(protein_acc) '
               )
        self.check_iso_and_peptide_relations(sql, proteincentric=True, nogroup=True)

    def check_iso_and_peptide_relations(self, sql, proteincentric=False, nogroup=False, flr=False):
        valsql = ('SELECT ps.sequence, bs.set_name, '
               'ppq.quant, pf.fdr, pp.pep, flr.flr '
               'FROM peptide_sequences AS ps '
               'JOIN biosets AS bs '
               'JOIN peptide_precur_quanted AS ppq USING(pep_id) '
               'LEFT OUTER JOIN ptm_flr AS flr USING(pep_id) '
               'JOIN peptide_fdr AS pf USING(pep_id) '
               'JOIN peptide_pep AS pp USING(pep_id) ')
        fields = ['MS1 area (highest of all PSMs)', 'q-value', 'PEP']
        if flr:
            fields += ['PTM FLR']
        self.check_build_values(valsql, fields, 'Peptide sequence')
        isosql = ('SELECT ps.sequence, bs.set_name, ch.channel_name, '
               'iq.quantvalue, iq.amount_psms, fqp.amount_psms '
               'FROM peptide_sequences AS ps '
               'JOIN biosets AS bs '
               'JOIN peptide_iso_quanted AS iq USING(pep_id) '
               'JOIN pepquant_channels AS ch USING(channel_id) '
               'JOIN peptide_iso_fullpsms AS fqp USING(pep_id) '
               )
        self.check_built_isobaric(isosql, 'Peptide sequence')
        expected, psm_id = {}, None
        for rec in self.get_values_from_db(self.dbfile, sql):
            rec = [rec[0]] + ['NA' if x is None else x for x in rec[1:]]
            try:
                expected[rec[0]]['psms'].add(rec[1])
            except KeyError:
                expected[rec[0]] = {'psms': set([rec[1]]),
                                    'pgroups': set([rec[2]]),
                                    'descriptions': set([rec[3]]),
                                    'genes': set([rec[4]]),
                                    'assoc': set([rec[5]]),
                                    'cover': set([rec[6]]),
                                    'pseqs': {rec[2]: rec[7]},
                                    }
            else:
                expected[rec[0]]['pgroups'].add(rec[2])
                expected[rec[0]]['descriptions'].add(rec[3])
                expected[rec[0]]['genes'].add(rec[4])
                expected[rec[0]]['assoc'].add(rec[5])
                expected[rec[0]]['cover'].add(rec[6])
                expected[rec[0]]['pseqs'][rec[2]] = rec[7]
        for line in self.tsv_generator(self.resultfn):
            pepseq = line['Peptide sequence']
            barepep = re.sub('[^A-Z]', '', pepseq)
            self.assertEqual(line['Bare peptide'], barepep)
            if proteincentric or nogroup:
                proteins = line['Protein(s)'].split(';')
                for pacc, ss in zip(proteins, line['Protein start/stop site(s)'].split(';')):
                    start, stop = [int(x) for x in ss.split('-')]
                    ix = expected[pepseq]['pseqs'][pacc].index(barepep)
                    self.assertEqual(start, ix + 1)
                    self.assertEqual(stop, ix + len(barepep))
                
            if proteincentric:
                self.assertEqual(set(line['Protein(s)'].split(';')),
                                 expected[pepseq]['pgroups'])
                if not nogroup:
                    self.assertEqual(
                        set(line['Description(s)'].split(';')),
                        set([y for x in
                             expected[pepseq]['descriptions']
                             for y in x.split(';')]))
                    rescovs = sorted(line['Coverage(s)'].split(';'))
                    expcovs = sorted(expected[pepseq]['cover'])
                    for rescov, expcov in zip(rescovs, expcovs):
                        self.assertAlmostEqual(float(rescov), expcov)
            if not nogroup:
                self.assertEqual(set(line['Gene ID(s)'].split(';')),
                                 expected[pepseq]['genes'])
                self.assertEqual(set(line['Gene name(s)'].split(';')),
                                 expected[pepseq]['assoc'])


class TestProteinMerge(basetests.MergeTest):

    def proteins(self, centric='proteincentric', cutoff=False):
        self.run_command(self.options)
        sql = ('SELECT p.protein_acc, bs.set_name, pf.fdr, pp.quant FROM proteins AS p '
               'JOIN biosets AS bs '
               'JOIN protein_fdr AS pf USING(pacc_id) '
               'JOIN protein_precur_quanted AS pp USING(pacc_id) '
               )
        self.check_build_values(sql, ['q-value', 'MS1 precursor area'],
                'Protein ID', cutoff)
        sql = ('SELECT p.protein_acc, bs.set_name, pc.channel_name, '
               'pi.quantvalue, pi.amount_psms, fqp.amount_psms, pf.fdr '
               'FROM proteins AS p '
               'JOIN biosets AS bs '
               'JOIN protein_iso_quanted AS pi USING(pacc_id) '
               'JOIN protein_iso_fullpsms AS fqp USING(pacc_id) '
               'JOIN protquant_channels AS pc USING(channel_id) '
               'JOIN protein_fdr AS pf WHERE pf.prottable_id=pc.prottable_id '
               'AND pf.pacc_id=p.pacc_id'
               )
        self.check_built_isobaric(sql, 'Protein ID', cutoff=cutoff)
        sql = """
        SELECT p.protein_acc,GROUP_CONCAT(g.gene_acc),GROUP_CONCAT(aid.assoc_id),pd.description,pcov.coverage
        FROM protein_group_master AS pgm 
        INNER JOIN proteins AS p ON p.pacc_id=pgm.pacc_id 
        LEFT OUTER JOIN genename_proteins AS gnp ON gnp.pacc_id=p.pacc_id
        LEFT OUTER JOIN associated_ids AS aid ON aid.gn_id=gnp.gn_id
        LEFT OUTER JOIN ensg_proteins AS egp ON egp.pacc_id=p.pacc_id
        LEFT OUTER JOIN genes AS g ON g.gene_id=egp.gene_id
        LEFT OUTER JOIN prot_desc AS pd USING(pacc_id) 
        JOIN protein_coverage AS pcov ON pcov.protein_acc=p.protein_acc
        GROUP BY pgm.pacc_id
        """
        psm_sql = """
        SELECT x.protein_acc, pp.psm_id, ps.sequence
        FROM protein_group_master AS pgm
        INNER JOIN proteins AS x ON x.pacc_id=pgm.pacc_id
        JOIN protein_psm AS pp ON pp.protein_acc=x.protein_acc
        JOIN psms USING(psm_id)
        JOIN mzml USING(spectra_id)
        JOIN mzmlfiles USING(mzmlfile_id)
        JOIN biosets AS bs USING(set_id)
        JOIN peptide_sequences AS ps USING(pep_id)
        WHERE bs.set_name='Set1'
        """
        self.check_protein_data(centric, sql, psm_sql)

    def test_proteins(self):
        self.infilename = 'proteins.txt'
        self.proteins()

    def test_mergecutoff(self):
        # FIXME need worse peptides, all FDR is 0.0, cutoff not tested!!
        # test is in check_built_isobaric
        self.infilename = 'proteins.txt'
        cutoff= 0.0001
        self.options.extend(['--mergecutoff', str(cutoff)])
        self.proteins(cutoff=cutoff)

    def test_nopsmnrs(self):
        """Given input tables with NO amount psm information, output a nice gene
        table without those columns"""
        isosql = """
        SELECT g.gene_acc, bs.set_name, pc.channel_name,
        pi.quantvalue, pi.amount_psms FROM genes AS g
        JOIN biosets AS bs
        JOIN gene_iso_quanted AS pi USING(gene_id)
        JOIN genequant_channels AS pc USING(channel_id)
        """
        self.ensgcentric_in_memory('ensg_nopsms.txt', isosql, nopsms=True)

    def test_nogroups(self):
        self.options.append('--no-group-annotation')
        self.infilename = 'proteins.txt'
        self.proteins(centric='protnogroup')

    def test_ensgcentric(self):
        isosql = """
        SELECT g.gene_acc, bs.set_name, pc.channel_name,
        pi.quantvalue, pi.amount_psms, fqp.amount_psms FROM genes AS g
        JOIN biosets AS bs
        JOIN gene_iso_quanted AS pi USING(gene_id)
        JOIN gene_iso_fullpsms AS fqp USING(gene_id)
        JOIN genequant_channels AS pc USING(channel_id)
        """
        self.ensgcentric_in_memory('ensg.txt', isosql)

    def ensgcentric_in_memory(self, infile, isosql, nopsms=False):
        self.infilename = infile
        self.run_command([*self.options, '--in-memory'])
        sql = ('SELECT g.gene_acc, bs.set_name, gf.fdr, gp.quant '
               'FROM genes AS g '
               'JOIN biosets AS bs '
               'JOIN gene_fdr AS gf USING(gene_id) '
               'JOIN gene_precur_quanted AS gp USING(gene_id) '
               )
        self.check_build_values(sql, ['q-value', 'MS1 precursor area'],
                'Gene ID')
        self.check_built_isobaric(isosql, 'Gene ID', nopsms)
        sql = """
        SELECT g.gene_acc, GROUP_CONCAT(aid.assoc_id),GROUP_CONCAT(p.protein_acc),pd.description
        FROM genes AS g
        INNER JOIN ensg_proteins AS egp ON egp.gene_id=g.gene_id
        LEFT OUTER JOIN genename_proteins AS gnp ON gnp.pacc_id=egp.pacc_id
        JOIN associated_ids AS aid ON aid.gn_id=gnp.gn_id
        JOIN proteins AS p ON p.pacc_id=gnp.pacc_id
        JOIN prot_desc AS pd ON pd.pacc_id=p.pacc_id
        GROUP BY g.gene_acc
        """
        psm_sql = """
        SELECT x.gene_acc, pp.psm_id, ps.sequence
        FROM genes AS x
        INNER JOIN ensg_proteins AS egp ON egp.gene_id=x.gene_id
        INNER JOIN proteins AS p ON p.pacc_id=egp.pacc_id
        JOIN protein_psm AS pp ON pp.protein_acc=p.protein_acc
        JOIN psms USING(psm_id)
        JOIN mzml USING(spectra_id)
        JOIN mzmlfiles USING(mzmlfile_id)
        JOIN biosets AS bs USING(set_id)
        JOIN peptide_sequences AS ps USING(pep_id)
        WHERE bs.set_name='Set1'
        """
        self.check_protein_data('genecentric', sql, psm_sql)

    def test_genenamecentric(self):
        self.infilename = 'genenames.txt'
        self.run_command(self.options)
        sql = ('SELECT ai.assoc_id, bs.set_name, gf.fdr, gp.quant '
               'FROM associated_ids AS ai '
               'JOIN biosets AS bs '
               'JOIN assoc_fdr AS gf USING(gn_id) '
               'JOIN assoc_precur_quanted AS gp USING(gn_id) '
               )
        self.check_build_values(sql, ['q-value', 'MS1 precursor area'], 
                'Gene Name')
        sql = ('SELECT ai.assoc_id, bs.set_name, pc.channel_name, '
               'pi.quantvalue, pi.amount_psms, fqp.amount_psms '
               'FROM associated_ids AS ai '
               'JOIN biosets AS bs '
               'JOIN assoc_iso_quanted AS pi USING(gn_id) '
               'JOIN assoc_iso_fullpsms AS fqp USING(gn_id) '
               'JOIN genequant_channels AS pc USING(channel_id) '
               )
        self.check_built_isobaric(sql, 'Gene Name')
        sql = """
        SELECT aid.assoc_id, GROUP_CONCAT(g.gene_acc), GROUP_CONCAT(p.protein_acc),pd.description
        FROM associated_ids AS aid
        INNER JOIN genename_proteins AS gnp ON gnp.gn_id=aid.gn_id
        LEFT OUTER JOIN ensg_proteins AS egp ON egp.pacc_id=gnp.pacc_id
        JOIN genes AS g ON g.gene_id=egp.gene_id
        JOIN proteins AS p ON p.pacc_id=gnp.pacc_id
        JOIN prot_desc AS pd ON pd.pacc_id=p.pacc_id
        GROUP BY aid.assoc_id
        """
        psm_sql = """
        SELECT x.assoc_id, pp.psm_id, ps.sequence
        FROM associated_ids AS x
        INNER JOIN genename_proteins AS gnp ON gnp.gn_id=x.gn_id
        INNER JOIN proteins AS p ON p.pacc_id=gnp.pacc_id
        JOIN protein_psm AS pp ON pp.protein_acc=p.protein_acc
        JOIN psms USING(psm_id)
        JOIN mzml USING(spectra_id)
        JOIN mzmlfiles USING(mzmlfile_id)
        JOIN biosets AS bs USING(set_id)
        JOIN peptide_sequences AS ps USING(pep_id)
        WHERE bs.set_name='Set1'
        """
        self.check_protein_data('assoccentric', sql, psm_sql) 

    def check_peptide_relations(self, sql):
        expected, psm_id = {}, None
        for rec in self.get_values_from_db(self.dbfile, sql):
            rec = [rec[0]] + ['NA' if x is None else x for x in rec[1:]]
            try:
                expected[rec[0]]['psms'].add(rec[1])
            except KeyError:
                expected[rec[0]] = {'psms': set([rec[1]]),
                                    'pgroups': set([rec[2]]),
                                    'descriptions': set([rec[3]]),
                                    'genes': set([rec[4]]),
                                    'assoc': set([rec[5]]),
                                    'cover': set([str(rec[6])]),
                                    }
            else:
                expected[rec[0]]['pgroups'].add(rec[2])
                expected[rec[0]]['descriptions'].add(rec[3])
                expected[rec[0]]['genes'].add(rec[4])
                expected[rec[0]]['assoc'].add(rec[5])
                expected[rec[0]]['cover'].add(str(rec[6]))
        for line in self.tsv_generator(self.resultfn):
            self.assertEqual(set(line['Protein(s)'].split(';')),
                             expected[line['Peptide sequence']]['pgroups'])
            self.assertEqual(
                set(line['Description(s)'].split(';')),
                set([y for x in
                     expected[line['Peptide sequence']]['descriptions']
                     for y in x.split(';')]))
            self.assertEqual(set(line['Gene ID(s)'].split(';')),
                             expected[line['Peptide sequence']]['genes'])
            self.assertEqual(set(line['Gene name(s)'].split(';')),
                             expected[line['Peptide sequence']]['assoc'])
            self.assertEqual(set([line['Coverage(s)']]),
                             expected[line['Peptide sequence']]['cover'])
