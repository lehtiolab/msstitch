from app.lookups.sqlite.base import ResultLookupInterface


MASTER_INDEX = 1
PROTEIN_ACC_INDEX = 2
PEPTIDE_COUNT_INDEX = 3
PSM_COUNT_INDEX = 4
PROTEIN_SCORE_INDEX = 5
EVIDENCE_LVL_INDEX = 6
COVERAGE_INDEX = 7


class ProteinGroupDB(ResultLookupInterface):
    def add_tables(self):
        self.create_tables(['proteins', 'protein_psm',
                            'protein_evidence', 'protein_seq',
                            'protein_coverage', 'protein_group_master',
                            'protein_group_content', 'psm_protein_groups',
                            'prot_desc'])

    def store_proteins(self, proteins, evidence_lvls, sequences=False):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO proteins(protein_acc) '
            'VALUES(?)', proteins)
        cursor.executemany(
            'INSERT INTO protein_evidence(protein_acc, evidence_lvl) '
            'VALUES(?, ?)', evidence_lvls)
        if sequences:
            cursor.executemany(
                'INSERT INTO protein_seq(protein_acc, sequence) '
                'VALUES(?, ?)', sequences)
        self.conn.commit()

    def store_descriptions(self, descriptions):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO prot_desc(protein_acc, description) '
            'VALUES(?, ?)', descriptions)
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

    def index_protein_peptides(self):
        self.index_column('protein_index', 'protein_psm', 'protein_acc')
        self.index_column('protpsmid_index', 'protein_psm', 'psm_id')
        self.index_column('protdesc_index', 'prot_desc', 'protein_acc')

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

    def update_master_proteins(self, new_masters):
        cur = self.get_cursor()
        sql = 'UPDATE protein_group_master SET protein_acc=? WHERE master_id=?'
        cur.executemany(sql, new_masters)
        self.conn.commit()

    def get_master_ids(self, invert=False):
        cur = self.get_cursor()
        cur.execute('SELECT protein_acc, master_id FROM protein_group_master')
        if invert:
            return {p_acc: master_id for (p_acc, master_id) in cur}
        else:
            return {master_id: p_acc for (p_acc, master_id) in cur}

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
        sql = ('SELECT p.protein_acc, ps.sequence, pp.psm_id, psms.sequence '
               'FROM proteins AS p '
               'JOIN protein_seq AS ps USING(protein_acc) '
               'JOIN protein_psm AS pp USING(protein_acc) '
               'JOIN psms AS psms USING(psm_id)'
               )
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_master_contentproteins_psms(self):
        sql = ('SELECT ppg.master_id, ppg.psm_id, pp.protein_acc, p.sequence, '
               'p.score '
               'FROM psm_protein_groups AS ppg '
               'JOIN protein_psm AS pp USING(psm_id) '
               'JOIN psms AS p USING(psm_id) '
               'ORDER BY ppg.master_id'
               )
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_all_master_psms(self):
        sql = ('SELECT master_id, psm_id '
               'FROM psm_protein_groups '
               'ORDER BY master_id'
               )
        cursor = self.get_cursor()
        return ((master, psm)
                for master, psm in cursor.execute(sql).fetchall())

    def get_all_psms_proteingroups(self, coverage):
        fields = ['pr.rownr', 'ppg.master_id', 'pgc.protein_acc',
                  'pgc.peptide_count', 'pgc.psm_count', 'pgc.protein_score',
                  'pev.evidence_lvl']
        joins = [('psm_protein_groups', 'ppg', 'psm_id'),
                 ('protein_group_content', 'pgc', 'master_id'),
                 ('protein_evidence', 'pev', 'protein_acc')]
        if coverage:
            fields.append('pc.coverage')
            joins.append(('protein_coverage', 'pc', 'protein_acc'))
        join_sql = '\n'.join(['JOIN {0} AS {1} USING({2})'.format(
            j[0], j[1], j[2]) for j in joins])
        sql = 'SELECT {0} FROM psmrows AS pr {1} ORDER BY pr.rownr'.format(
            ', '.join(fields), join_sql)
        cursor = self.get_cursor()
        return cursor.execute(sql)
