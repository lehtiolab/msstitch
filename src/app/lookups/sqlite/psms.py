from app.lookups.sqlite.base import ResultLookupInterface


class PSMDB(ResultLookupInterface):
    def add_tables(self, tabletypes):
        self.create_tables(['psms', 'psmrows', 'peptide_sequences',
                            'proteins', 'protein_evidence', 'protein_seq',
                            'prot_desc', 'protein_psm', 'genes',
                            'associated_ids'])
        if 'proteingroup' in tabletypes:
            self.create_tables(['protein_coverage', 'protein_group_master',
                                'protein_group_content', 'psm_protein_groups'])

    def store_fasta(self, prot, evids, seq, desc, ensg, symbols):
        self.store_proteins(prot, evidence_lvls=evids, sequences=seq)
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO prot_desc(protein_acc, description) '
            'VALUES(?, ?)', desc)
        self.conn.commit()
        self.index_column('protdesc_index', 'prot_desc', 'protein_acc')
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO genes(gene_acc, protein_acc) VALUES(?, ?)', ensg)
        self.conn.commit()
        self.index_column('gene_index', 'genes', 'protein_acc')
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO associated_ids(assoc_id, protein_acc) VALUES(?, ?)',
            symbols)
        self.conn.commit()
        self.index_column('associd_index', 'associated_ids', 'protein_acc')

    def store_proteins(self, proteins, evidence_lvls=False, sequences=False):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO proteins(protein_acc) '
            'VALUES(?)', proteins)
        self.conn.commit()
        cursor = self.get_cursor()
        if evidence_lvls:
            cursor.executemany(
                'INSERT INTO protein_evidence(protein_acc, evidence_lvl) '
                'VALUES(?, ?)', evidence_lvls)
        if sequences:
            cursor.executemany(
                'INSERT INTO protein_seq(protein_acc, sequence) '
                'VALUES(?, ?)', sequences)
        self.conn.commit()
        self.index_column('proteins_index', 'proteins', 'protein_acc')
        self.index_column('evidence_index', 'protein_evidence', 'protein_acc')

    def get_protein_gene_map(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT p.protein_acc, g.gene_acc, aid.assoc_id, d.description '
            'FROM proteins AS p '
            'LEFT OUTER JOIN genes AS g ON p.protein_acc=g.protein_acc '
            'LEFT OUTER JOIN associated_ids AS aid '
            'ON p.protein_acc=aid.protein_acc '
            'LEFT OUTER JOIN prot_desc AS d ON p.protein_acc=d.protein_acc'
        )
        gpmap = {p_acc: {'gene': gene, 'symbol': sym, 'desc': desc}
                 for p_acc, gene, sym, desc in cursor}
        return gpmap

    def store_pepseqs(self, sequences):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO peptide_sequences(sequence) VALUES(?)', sequences)
        self.conn.commit()

    def store_psms(self, psms):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO psms(psm_id, pep_id, score, spectra_id) '
            'VALUES(?, ?, ?, ?)', ((psm['psm_id'], psm['seq'], psm['score'],
                                    psm['spec_id'])
                                   for psm in psms))
        cursor.executemany(
            'INSERT INTO psmrows(psm_id, rownr) VALUES(?, ?)',
            ((psm['psm_id'], psm['rownr']) for psm in psms))
        self.conn.commit()

    def store_peptides_proteins(self, allpepprot, psmids_to_store):
        ppmap = {psm_id: allpepprot[psm_id] for psm_id in psmids_to_store}
        prot_psm_ids = ((prot_acc, psm_id)
                        for psm_id, prots in ppmap.items()
                        for prot_acc in prots)
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO protein_psm(protein_acc, psm_id)'
            ' VALUES (?, ?)', prot_psm_ids)
        self.conn.commit()

    def get_peptide_seq_map(self):
        cursor = self.get_cursor()
        seqs = cursor.execute('SELECT pep_id, sequence FROM peptide_sequences')
        return {seq: pepid for pepid, seq in seqs.fetchall()}

    def index_psms(self):
        self.index_column('psmid_index', 'psms', 'psm_id')
        self.index_column('psmspecid_index', 'psms', 'spectra_id')
        self.index_column('psmrowid_index', 'psmrows', 'psm_id')
        self.index_column('psmrow_index', 'psmrows', 'rownr')
        self.index_column('pepseq_index', 'peptide_sequences', 'sequence')
        self.index_column('pepid_index', 'peptide_sequences', 'pep_id')
        self.index_column('psmspepid_index', 'psms', 'pep_id')

    def index_protein_peptides(self):
        self.index_column('protein_index', 'protein_psm', 'protein_acc')
        self.index_column('protpsmid_index', 'protein_psm', 'psm_id')

    def get_exp_spectra_data_rows(self):
        cursor = self.get_cursor()
        return cursor.execute('SELECT pr.rownr, bs.set_name, sp.retention_time, '
                              'iit.ion_injection_time, im.ion_mobility '
                              'FROM psmrows AS pr '
                              'JOIN psms AS p USING(psm_id) '
                              'JOIN mzml AS sp USING(spectra_id) '
                              'LEFT OUTER JOIN ioninjtime AS iit USING(spectra_id) '
                              'LEFT OUTER JOIN ionmob AS im USING(spectra_id) '
                              'JOIN mzmlfiles as mf USING(mzmlfile_id) '
                              'JOIN biosets AS bs USING(set_id) '
                              'ORDER BY pr.rownr')

    def store_masters(self, allmasters, psm_masters):
        allmasters = ((x,) for x in allmasters)
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO protein_group_master(protein_acc) VALUES(?)',
            allmasters)
        master_ids = self.get_master_ids()
        psms = ((psm_id, master_ids[master])
                for psm_id, masters in psm_masters.items()
                for master in masters)
        cursor.executemany(
            'INSERT INTO psm_protein_groups(psm_id, master_id) '
            'VALUES(?, ?)', psms)
        self.conn.commit()
        self.index_column('psm_pg_index', 'psm_protein_groups', 'master_id')
        self.index_column('psm_pg_psmid_index', 'psm_protein_groups', 'psm_id')

    def update_master_proteins(self, new_masters):
        cur = self.get_cursor()
        sql = 'UPDATE protein_group_master SET protein_acc=? WHERE master_id=?'
        cur.executemany(sql, new_masters)
        self.conn.commit()

    def get_master_ids(self):
        cur = self.get_cursor()
        cur.execute('SELECT protein_acc, master_id FROM protein_group_master')
        return {p_acc: master_id for (p_acc, master_id) in cur}

    def store_coverage(self, coverage):
        cursor = self.get_cursor()
        sql = ('INSERT INTO protein_coverage(protein_acc, coverage) '
               'VALUES(?, ?)')
        cursor.executemany(sql, coverage)
        self.conn.commit()
        self.index_column('cov_index', 'protein_coverage', 'protein_acc')

    def store_protein_group_content(self, protein_groups):
        cursor = self.get_cursor()
        cursor.executemany('INSERT INTO protein_group_content('
                           'protein_acc, master_id, peptide_count, '
                           'psm_count, protein_score) '
                           'VALUES(?, ?, ?, ?, ?)', protein_groups)
        self.conn.commit()

    def index_protein_group_content(self):
        self.index_column('pgc_master_index', 'protein_group_content',
                          'master_id')

    def get_all_psm_protein_relations(self):
        sql = 'SELECT psm_id, protein_acc FROM protein_psm'
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_allpsms_masters(self):
        sql = ('SELECT pgm.protein_acc, pp.psm_id FROM protein_group_master '
               'AS pgm JOIN protein_psm as pp USING(protein_acc)')
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_proteins_for_peptide(self, psm_id):
        """Returns list of proteins for a passed psm_id"""
        protsql = self.get_sql_select(['protein_acc'], 'protein_psm')
        protsql = '{0} WHERE psm_id=?'.format(protsql)
        cursor = self.get_cursor()
        proteins = cursor.execute(protsql, psm_id).fetchall()
        return [x[0] for x in proteins]

    def get_protpepmap_from_proteins(self, proteins):
        pepsql = self.get_sql_select(['protein_acc', 'psm_id'],
                                     'protein_psm',
                                     distinct=True)
        pepsql = '{0} WHERE protein_acc {1}'.format(
            pepsql, self.get_inclause(proteins))
        cursor = self.get_cursor()
        protpeps = cursor.execute(pepsql, proteins).fetchall()
        outmap = {}
        for protein, peptide in protpeps:
            try:
                outmap[protein].append(peptide)
            except KeyError:
                outmap[protein] = [peptide]
        return outmap

    def get_all_proteins_psms_seq(self):
        sql = ('SELECT p.protein_acc, ps.sequence, pp.psm_id, peps.sequence '
               'FROM proteins AS p '
               'JOIN protein_seq AS ps USING(protein_acc) '
               'JOIN protein_psm AS pp USING(protein_acc) '
               'JOIN psms AS psms USING(psm_id) '
               'JOIN peptide_sequences AS peps USING(pep_id)'
               )
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_protein_psm_records(self):
        sql = ('SELECT protein_acc, psm_id FROM protein_psm')
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_protein_group_candidates(self):
        sql = ('SELECT pgm.master_id, pgm.psm_id, pp.protein_acc, '
               'peps.sequence, p.score, pev.evidence_lvl, pc.coverage '
               'FROM psm_protein_groups AS pgm '
               'JOIN protein_psm AS pp USING(psm_id) '
               'JOIN psms AS p USING(psm_id) '
               'JOIN peptide_sequences AS peps USING(pep_id) '
               'LEFT OUTER JOIN protein_evidence AS pev '
               'ON pev.protein_acc=pp.protein_acc '
               'LEFT OUTER JOIN protein_coverage AS pc '
               'ON pc.protein_acc=pp.protein_acc '
               'ORDER BY pgm.master_id'
               )
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def check_evidence_tables(self):
        """Returns True if there are records in evidence
        tables, otherwise returns False"""
        return self.check_table('protein_evidence') is not False

    def check_table(self, tablename):
        cursor = self.get_cursor()
        cursor.execute('SELECT * FROM {} LIMIT 10'.format(tablename))
        if len(cursor.fetchall()) > 0:
            return True
        return False

    def get_all_psms_proteingroups(self, evidence):
        fields = ['pr.rownr', 'pgm.protein_acc', 'pgc.protein_acc',
                  'pgc.peptide_count', 'pgc.psm_count', 'pgc.protein_score',
                  'pc.coverage']
        joins = [('psm_protein_groups', 'ppg', 'psm_id'),
                 ('protein_group_master', 'pgm', 'master_id'),
                 ('protein_group_content', 'pgc', 'master_id'),
                 ]
        join_sql = '\n'.join(['JOIN {0} AS {1} USING({2})'.format(
            j[0], j[1], j[2]) for j in joins])
        specialjoin = []
        if evidence:
            specialjoin = [('protein_evidence', 'pev')]
            fields.append('pev.evidence_lvl')
        specialjoin.append(('protein_coverage', 'pc'))
        specialjoin = '\n'.join(['JOIN {0} AS {1} ON '
                                 'pgc.protein_acc={1}.protein_acc'.format(
                                     j[0], j[1]) for j in specialjoin])
        join_sql = '{} {}'.format(join_sql, specialjoin)
        sql = 'SELECT {0} FROM psmrows AS pr {1} ORDER BY pr.rownr'.format(
            ', '.join(fields), join_sql)
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def select_all_psm_quants(self, isobaric=False, precursor=False):
        selects = ['pr.rownr']
        joins = ['JOIN psms USING(psm_id)', 'JOIN mzml USING(spectra_id)']
        sqlfields, fieldcount = {}, 1
        if isobaric:
            selects.extend(['ic.channel_name', 'iq.intensity'])
            joins.extend(['JOIN isobaric_quant AS iq USING(spectra_id)',
                          'JOIN isobaric_channels AS ic USING(channel_id)'])
            sqlfields['isochan'] = fieldcount
            sqlfields['isoquant'] = fieldcount + 1
            fieldcount += 2
        if precursor:
            selects.extend(['pq.intensity'])
            joins.extend(['LEFT OUTER JOIN ms1_align USING(spectra_id)',
                          'LEFT OUTER JOIN ms1_quant AS pq USING(feature_id)'])
            sqlfields['precursor'] = fieldcount
        if not precursor and not isobaric:
            raise RuntimeError('Cannot add quantification data, neither '
                               'isobaric, nor precursor have been specified.')
        sql = ('SELECT {} FROM psmrows as pr {} '
               'ORDER BY pr.rownr'.format(', '.join(selects), ' '.join(joins)))
        cursor = self.get_cursor()
        return cursor.execute(sql), sqlfields

    def get_all_quantmaps(self):
        """Returns all unique quant channels from lookup as list"""
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT channel_name FROM isobaric_channels')
        return cursor.fetchall()
