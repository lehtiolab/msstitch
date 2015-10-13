from app.lookups.sqlite.base import ResultLookupInterface


class ProtPepTable(ResultLookupInterface):
    table_map = {'protein': {'fntable': 'protein_tables',
                             'feattable': 'proteins',
                             'prectable': 'protein_precur_quanted',
                             'fdrtable': 'protein_fdr',
                             'peptable': 'protein_pep',
                             },
                 'peptide': {'fntable': 'peptide_tables',
                             'feattable': 'peptide_sequences',
                             'prectable': 'peptide_precur_quanted',
                             'fdrtable': 'peptide_fdr',
                             'peptable': 'peptide_pep',
                             }
                 }

    def get_all_poolnames(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT set_name, set_id FROM biosets')
        return cursor

    def store_table_files(self, tables):
        self.store_many(
            'INSERT INTO {}(set_id, filename) VALUES(?, ?)'.format(
                self.table_map[self.datatype]['fntable']),
            tables)

    def get_tablefn_map(self):
        table = self.table_map[self.datatype]['fntable']
        cursor = self.get_cursor()
        cursor.execute('SELECT * FROM {} '.format(table))
        return {fn: table_id for (table_id, setid, fn) in cursor}

    def get_feature_map(self):
        columns = {'protein': ['pacc_id', 'protein_acc'],
                   'peptide': ['pep_id', 'sequence']
                   }
        table = self.table_map[self.datatype]['feattable']
        columns = columns[self.datatype]
        cursor = self.get_cursor()
        cursor.execute('SELECT {}, {} FROM {}'.format(columns[0], columns[1],
                                                      table))
        return {acc: table_id for (table_id, acc) in cursor}

    def store_singlecol(self, tablekey, vals):
        table = self.table_map[self.datatype][tablekey]
        cols = self.colmap[table]
        self.store_many('INSERT INTO {}({}, {}, {}) '
                        'VALUES (?, ?, ?)'.format(table, cols[0], cols[1],
                                                  cols[2]), vals)

    def store_precursor_quants(self, quants):
        self.store_singlecol('prectable', quants)

    def store_fdr(self, fdr):
        self.store_singlecol('fdrtable', fdr)

    def store_pep(self, pep):
        self.store_singlecol('peptable', pep)

    def update_selects(self, selectmap, fields, fieldcount):
        selectmap.update({field: i + fieldcount
                          for i, field in enumerate(fields)})
        fieldcount = max(selectmap.values()) + 1
        return selectmap, fieldcount

    def get_proteins_psms_genecentric(self):
        fields = ['p.protein_acc', 'sets.set_name',
                  'pep.sequence', 'psm.psm_id']
        firstjoin = ('protein_psm', 'pp', 'protein_acc')
        return self.get_proteins_psms('protein', fields, firstjoin)

    def get_proteins_psms_pgrouped(self):
        fields = ['p.protein_acc', 'sets.set_name',
                  'pep.sequence', 'psm.psm_id', 'pd.description',
                  'pcov.coverage']
        leftjoins = ('LEFT OUTER JOIN prot_desc AS pd USING(protein_acc) '
                     'LEFT OUTER JOIN protein_coverage '
                     'AS pcov USING(protein_acc)'
                     )
        firstjoin = ('psm_protein_groups', 'ppg', 'master_id')
        return self.get_proteins_psms('protein_group_master', fields,
                                      firstjoin, leftjoins)

    def get_proteins_psms(self, firsttable, fields, firstjoin,
                          leftjoins=False):
        joins = [firstjoin]
        joins.extend([('psm_protein_groups', 'ppg', 'master_id'),
                      ('psms', 'psm', 'psm_id'),
                      ('peptide_sequences', 'pep', 'pep_id'),
                      ('mzml', 'sp', 'spectra_id'),
                      ('mzmlfiles', 'mzfn', 'mzmlfile_id'),
                      ('biosets', 'sets', 'set_id'),
                      ])
        join_sql = ' '.join(['JOIN {} AS {} USING({})'.format(
            j[0], j[1], j[2]) for j in joins])
        if leftjoins:
            join_sql = '{} {}'.format(join_sql, leftjoins)
        sql = ('SELECT {} FROM {}'
               'AS p'.format(', '.join(fields), firsttable))
        sql = '{} {} ORDER BY p.protein_acc, sets.set_name'.format(sql,
                                                                   join_sql)
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_sql_joins_mergetable(self, sql, joins, pep_or_prot):
        protein_j_cols = {'p': 'pacc_id', 'pt': 'prottable_id'}
        peptide_j_cols = {'p': 'pep_id', 'pt': 'peptable_id'}
        colpick = {'peptide': peptide_j_cols, 'protein': protein_j_cols}
        join_cols = {'pc': 'channel_id'}
        join_cols.update(colpick[pep_or_prot])
        if joins:
            joinsql = ''
            for j in joins:
                joincmd = 'JOIN'
                if j[-1] is True:
                    joincmd = 'LEFT OUTER {}'.format(joincmd)
                joinmatchsql = ' AND '.join(['{0}.{1}={2}.{1}'.format(
                    j[1], join_cols[jj], jj) for jj in j[2]])
                joinsql = '{4} {0} {1} AS {2} ON {3}'.format(
                    joincmd, j[0], j[1], joinmatchsql, joinsql)
            sql = '{} {}'.format(sql, joinsql)
        return sql

    def get_merged_features(self, sql):
        return self.execute_sql(sql)
