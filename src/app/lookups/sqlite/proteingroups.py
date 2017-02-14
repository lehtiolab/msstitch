from app.lookups.sqlite.base import ResultLookupInterface


# Indices that belong to positions of these features in output from
# function get_all_psms_proteingroups:
MASTER_INDEX = 1
PROTEIN_ACC_INDEX = 2
PEPTIDE_COUNT_INDEX = 3
PSM_COUNT_INDEX = 4
PROTEIN_SCORE_INDEX = 5
COVERAGE_INDEX = 6
EVIDENCE_LVL_INDEX = 7


class ProteinGroupDB(ResultLookupInterface):
    def add_tables(self):
        self.create_tables(['protein_coverage', 'protein_group_master',
                            'protein_group_content', 'psm_protein_groups'])

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
