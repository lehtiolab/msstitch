import sqlite3


mslookup_tables = {'biosets': ['set_id INTEGER PRIMARY KEY',
                               'set_name TEXT'],
                   'mzmlfiles': ['mzmlfile_id INTEGER PRIMARY KEY',
                                 'mzmlfilename TEXT',
                                 'set_id INTEGER',
                                 'FOREIGN KEY(set_id)'
                                 'REFERENCES biosets'
                                 ],
                   'mzml': ['spectra_id TEXT PRIMARY KEY',
                            'mzmlfile_id INTEGER',
                            'scan_nr TEXT',
                            'charge INTEGER',
                            'mz DOUBLE',
                            'retention_time DOUBLE',
                            'ion_injection_time DOUBLE',
                            'FOREIGN KEY(mzmlfile_id)'
                            'REFERENCES mzmlfiles'],
                   'isobaric_channels': ['channel_id INTEGER PRIMARY KEY',
                                         'channel_name TEXT'],
                   'isobaric_quant': ['spectra_id TEXT',
                                      'channel_id INTEGER',
                                      'intensity REAL',
                                      'FOREIGN KEY(spectra_id)'
                                      'REFERENCES mzml',
                                      'FOREIGN KEY(channel_id)'
                                      'REFERENCES isobaric_channels'
                                      ],
                   # ms1_quant has no spectra_id reference since it contains
                   # features and Im not sure if they can be linked to
                   # spectra_ids or that retention time is averaged
                   'ms1_quant': ['feature_id INTEGER PRIMARY KEY',
                                 'mzmlfile_id INTEGER',
                                 'retention_time REAL', 'mz REAL',
                                 'charge INTEGER', 'intensity REAL',
                                 'FOREIGN KEY(mzmlfile_id)'
                                 'REFERENCES mzmlfiles'],
                   'ms1_align': ['spectra_id TEXT',
                                 'feature_id INTEGER',
                                 'FOREIGN KEY(spectra_id) '
                                 'REFERENCES mzml '
                                 'FOREIGN KEY(feature_id) '
                                 'REFERENCES ms1_quant',
                                 ],
                   'peptide_sequences': ['pep_id INTEGER PRIMARY KEY',
                                         'sequence TEXT',
                                         ],
                   'psms': ['psm_id TEXT PRIMARY KEY NOT NULL',
                            'pep_id INTEGER',
                            'score TEXT',
                            'spectra_id TEXT',
                            'FOREIGN KEY(pep_id)'
                            'REFERENCES peptide_sequences '
                            'FOREIGN KEY(spectra_id)'
                            'REFERENCES mzml'
                            ],
                   'psmrows': ['psm_id TEXT',
                               'rownr INTEGER',
                               'FOREIGN KEY(psm_id) '
                               'REFERENCES psms(psm_id)'],
                   'proteins': ['pacc_id INTEGER PRIMARY KEY',
                                'protein_acc TEXT UNIQUE'],
                   'gene_tables': ['genetable_id INTEGER PRIMARY KEY',
                                   'set_id INTEGER',
                                   'filename TEXT',
                                   'FOREIGN KEY(set_id)'
                                   'REFERENCES biosets'],
                   'protein_tables': ['prottable_id INTEGER PRIMARY KEY',
                                      'set_id INTEGER',
                                      'filename TEXT',
                                      'FOREIGN KEY(set_id)'
                                      'REFERENCES biosets'],
                   'peptide_tables': ['peptable_id INTEGER PRIMARY KEY',
                                      'set_id INTEGER',
                                      'filename TEXT',
                                      'FOREIGN KEY(set_id)'
                                      'REFERENCES biosets'],
                   'pepquant_channels': ['channel_id INTEGER PRIMARY KEY',
                                         'peptable_id INTEGER',
                                         'channel_name TEXT',
                                         'amount_psms_name TEXT',
                                         'FOREIGN KEY(peptable_id) '
                                         'REFERENCES '
                                         'peptide_tables(peptable_id)'
                                         ],
                   'peptide_iso_quanted': ['peptidequant_id'
                                           'INTEGER PRIMARY KEY',
                                           'pep_id INTEGER',
                                           'channel_id INTEGER',
                                           'quantvalue REAL',
                                           'amount_psms INTEGER',
                                           'FOREIGN KEY(pep_id) '
                                           'REFERENCES '
                                           'peptide_sequences(pep_id) '
                                           'FOREIGN KEY(channel_id) '
                                           'REFERENCES '
                                           'pepquant_channels(channel_id)'
                                           ],
                   'peptide_precur_quanted':
                   ['pep_precquant_id INTEGER PRIMARY KEY',
                    'pep_id INTEGER',
                    'peptable_id INTEGER',
                    'quant REAL',
                    'FOREIGN KEY(pep_id) '
                    'REFERENCES peptide_sequences(pep_id) '
                    'FOREIGN KEY(peptable_id) '
                    'REFERENCES peptide_tables(peptable_id)'
                    ],
                   'peptide_fdr': ['pep_id INTEGER',
                                   'peptable_id INTEGER',
                                   'fdr DOUBLE',
                                   'FOREIGN KEY(pep_id) '
                                   'REFERENCES peptide_sequences(pep_id) '
                                   'FOREIGN KEY(peptable_id) '
                                   'REFERENCES '
                                   'peptide_tables(peptable_id)'
                                   ],
                   'peptide_pep': ['pep_id INTEGER',
                                   'peptable_id INTEGER',
                                   'pep DOUBLE',
                                   'FOREIGN KEY(pep_id) '
                                   'REFERENCES peptide_sequences(pep_id) '
                                   'FOREIGN KEY(peptable_id) '
                                   'REFERENCES '
                                   'peptide_tables(peptable_id)'
                                   ],
                   'protein_precur_quanted':
                   ['prot_precquant_id INTEGER PRIMARY KEY',
                    'pacc_id INTEGER',
                    'prottable_id INTEGER',
                    'quant REAL',
                    'FOREIGN KEY(pacc_id) '
                    'REFERENCES proteins(pacc_id) '
                    'FOREIGN KEY(prottable_id) '
                    'REFERENCES protein_tables(prottable_id)'
                    ],
                   'gene_precur_quanted':
                   ['gene_precquant_id INTEGER PRIMARY KEY',
                    'gene_id INTEGER',
                    'genetable_id INTEGER',
                    'quant REAL',
                    'FOREIGN KEY(gene_id) '
                    'REFERENCES genes(gene_id) '
                    'FOREIGN KEY(genetable_id) '
                    'REFERENCES gene_tables(genetable_id)'
                    ],
                   'assoc_precur_quanted':
                   ['gene_precquant_id INTEGER PRIMARY KEY',
                    'gene_id INTEGER',
                    'genetable_id INTEGER',
                    'quant REAL',
                    'FOREIGN KEY(gene_id) '
                    'REFERENCES associated_ids(gene_id) '
                    'FOREIGN KEY(genetable_id) '
                    'REFERENCES gene_tables(genetable_id)'
                    ],
                   'protein_iso_quanted': ['proteinquant_id '
                                           'INTEGER PRIMARY KEY',
                                           'pacc_id INTEGER',
                                           'channel_id INTEGER',
                                           'quantvalue REAL',
                                           'amount_psms INTEGER',
                                           'FOREIGN KEY(pacc_id) '
                                           'REFERENCES proteins(pacc_id) '
                                           'FOREIGN KEY(channel_id) '
                                           'REFERENCES '
                                           'protquant_channels(channel_id)'
                                           ],
                   'gene_iso_quanted': ['genequant_id '
                                        'INTEGER PRIMARY KEY',
                                        'gene_id INTEGER',
                                        'channel_id INTEGER',
                                        'quantvalue REAL',
                                        'amount_psms INTEGER',
                                        'FOREIGN KEY(gene_id) '
                                        'REFERENCES genes(gene_id) '
                                        'FOREIGN KEY(channel_id) '
                                        'REFERENCES '
                                        'genequant_channels(channel_id)'
                                        ],
                   'assoc_iso_quanted': ['genequant_id '
                                         'INTEGER PRIMARY KEY',
                                         'gene_id INTEGER',
                                         'channel_id INTEGER',
                                         'quantvalue REAL',
                                         'amount_psms INTEGER',
                                         'FOREIGN KEY(gene_id) '
                                         'REFERENCES associated_ids(gene_id) '
                                         'FOREIGN KEY(channel_id) '
                                         'REFERENCES '
                                         'genequant_channels(channel_id)'
                                         ],
                   'genequant_channels': ['channel_id INTEGER PRIMARY KEY',
                                          'genetable_id INTEGER',
                                          'channel_name TEXT',
                                          'amount_psms_name TEXT',
                                          'FOREIGN KEY(genetable_id) '
                                          'REFERENCES '
                                          'gene_tables(genetable_id)'
                                          ],
                   'protquant_channels': ['channel_id INTEGER PRIMARY KEY',
                                          'prottable_id INTEGER',
                                          'channel_name TEXT',
                                          'amount_psms_name TEXT',
                                          'FOREIGN KEY(prottable_id) '
                                          'REFERENCES '
                                          'protein_tables(prottable_id)'
                                          ],
                   'gene_probability': ['gene_id INTEGER',
                                        'genetable_id INTEGER',
                                        'probability DOUBLE',
                                        'FOREIGN KEY(gene_id) '
                                        'REFERENCES genes(gene_id) '
                                        'FOREIGN KEY(genetable_id) '
                                        'REFERENCES '
                                        'gene_tables(genetable_id)'
                                        ],
                   'assoc_probability': ['gene_id INTEGER',
                                         'genetable_id INTEGER',
                                         'probability DOUBLE',
                                         'FOREIGN KEY(gene_id) '
                                         'REFERENCES associated_ids(gene_id) '
                                         'FOREIGN KEY(genetable_id) '
                                         'REFERENCES '
                                         'gene_tables(genetable_id)'
                                         ],
                   'protein_probability': ['pacc_id INTEGER',
                                           'prottable_id INTEGER',
                                           'probability DOUBLE',
                                           'FOREIGN KEY(pacc_id) '
                                           'REFERENCES proteins(pacc_id) '
                                           'FOREIGN KEY(prottable_id) '
                                           'REFERENCES '
                                           'protein_tables(prottable_id)'
                                           ],
                   'gene_fdr': ['gene_id INTEGER',
                                'genetable_id INTEGER',
                                'fdr DOUBLE',
                                'FOREIGN KEY(gene_id) '
                                'REFERENCES genes(gene_id) '
                                'FOREIGN KEY(genetable_id) '
                                'REFERENCES '
                                'gene_tables(genetable_id)'
                                ],
                   'assoc_fdr': ['gene_id INTEGER',
                                 'genetable_id INTEGER',
                                 'fdr DOUBLE',
                                 'FOREIGN KEY(gene_id) '
                                 'REFERENCES associated_ids(gene_id) '
                                 'FOREIGN KEY(genetable_id) '
                                 'REFERENCES '
                                 'gene_tables(genetable_id)'
                                 ],
                   'protein_fdr': ['pacc_id INTEGER',
                                   'prottable_id INTEGER',
                                   'fdr DOUBLE',
                                   'FOREIGN KEY(pacc_id) '
                                   'REFERENCES proteins(pacc_id) '
                                   'FOREIGN KEY(prottable_id) '
                                   'REFERENCES '
                                   'protein_tables(prottable_id)'
                                   ],
                   'gene_pep': ['gene_id INTEGER',
                                'genetable_id INTEGER',
                                'pep DOUBLE',
                                'FOREIGN KEY(gene_id) '
                                'REFERENCES genes(gene_id) '
                                'FOREIGN KEY(genetable_id) '
                                'REFERENCES '
                                'gene_tables(genetable_id)'
                                ],
                   'assoc_pep': ['gene_id INTEGER',
                                 'genetable_id INTEGER',
                                 'pep DOUBLE',
                                 'FOREIGN KEY(gene_id) '
                                 'REFERENCES associated_ids(gene_id) '
                                 'FOREIGN KEY(genetable_id) '
                                 'REFERENCES '
                                 'gene_tables(genetable_id)'
                                 ],
                   'protein_pep': ['pacc_id INTEGER',
                                   'prottable_id INTEGER',
                                   'pep DOUBLE',
                                   'FOREIGN KEY(pacc_id) '
                                   'REFERENCES proteins(pacc_id) '
                                   'FOREIGN KEY(prottable_id) '
                                   'REFERENCES '
                                   'protein_tables(prottable_id)'
                                   ],
                   'protein_psm': ['protein_acc TEXT',
                                   'psm_id TEXT',
                                   'FOREIGN KEY(protein_acc) '
                                   'REFERENCES proteins(protein_acc) '
                                   'FOREIGN KEY(psm_id) '
                                   'REFERENCES psms(psm_id)'],
                   'protein_evidence': ['protein_acc TEXT',
                                        'evidence_lvl REAL',
                                        'FOREIGN KEY(protein_acc) '
                                        'REFERENCES '
                                        'proteins(protein_acc)'],
                   'protein_seq': ['protein_acc TEXT',
                                   'sequence TEXT',
                                   'FOREIGN KEY(protein_acc) '
                                   'REFERENCES '
                                   'proteins(protein_acc)'],
                   'protein_coverage': ['protein_acc TEXT',
                                        'coverage REAL',
                                        'FOREIGN KEY(protein_acc) '
                                        'REFERENCES '
                                        'proteins(protein_acc)'],
                   'protein_group_master': ['master_id INTEGER PRIMARY KEY',
                                            'protein_acc TEXT',
                                            'FOREIGN KEY(protein_acc) '
                                            'REFERENCES '
                                            'proteins(protein_acc)'],
                   'protein_group_content': ['protein_acc TEXT',
                                             'master_id INTEGER',
                                             'peptide_count INTEGER',
                                             'psm_count INTEGER',
                                             'protein_score INTEGER',
                                             'FOREIGN KEY(protein_acc) '
                                             'REFERENCES '
                                             'proteins(protein_acc) '
                                             'FOREIGN KEY(master_id) '
                                             'REFERENCES '
                                             'protein_group_master'
                                             '(master_id)'
                                             ],
                   'psm_protein_groups': ['psm_id TEXT',
                                          'master_id INTEGER',
                                          'FOREIGN KEY(psm_id) REFERENCES'
                                          ' psms(psm_id)',
                                          'FOREIGN KEY(master_id) REFERENCES'
                                          ' protein_group_master(master_id)'],
                   'genes': ['gene_id INTEGER PRIMARY KEY',
                             'gene_acc TEXT',
                             'protein_acc TEXT',
                             'FOREIGN KEY(protein_acc) '
                             'REFERENCES proteins(protein_acc)'],
                   'associated_ids': ['gene_id INTEGER PRIMARY KEY',
                                      'assoc_id TEXT',
                                      'protein_acc TEXT',
                                      'FOREIGN KEY(protein_acc) '
                                      'REFERENCES proteins(protein_acc)'],
                   'prot_desc': ['protein_acc TEXT',
                                 'description TEXT',
                                 'FOREIGN KEY(protein_acc) '
                                 'REFERENCES proteins(protein_acc)'],
                   'known_searchspace': ['seqs TEXT'],
                   'protein_peptides': ['seq TEXT', 'protid TEXT',
                                        'pos INTEGER'],
                   }


class DatabaseConnection(object):
    def __init__(self, fn=None):
        """SQLite connecting when given filename"""
        self.fn = fn
        if self.fn is not None:
            self.connect(self.fn)

    def get_fn(self):
        """Returns lookup filename"""
        return self.fn

    def create_tables(self, tables):
        """Creates database tables in sqlite lookup db"""
        cursor = self.get_cursor()
        for table in tables:
            columns = mslookup_tables[table]
            try:
                cursor.execute('CREATE TABLE {0}({1})'.format(
                    table, ', '.join(columns)))
            except sqlite3.OperationalError as error:
                print(error)
                print('Warning: Table {} already exists in database, will '
                      'add to existing tables instead of creating '
                      'new.'.format(table))
            else:
                self.conn.commit()

    def connect(self, fn):
        """SQLite connect method initialize db"""
        self.conn = sqlite3.connect(fn)
        cur = self.get_cursor()
        cur.execute('PRAGMA page_size=4096')
        cur.execute('PRAGMA FOREIGN_KEYS=ON')
        cur.execute('PRAGMA cache_size=10000')
        cur.execute('PRAGMA journal_mode=MEMORY')

    def get_cursor(self):
        """Quickly get cursor, abstracting connection"""
        return self.conn.cursor()

    def close_connection(self):
        """Close connection to db, abstracts connection object"""
        self.conn.close()

    def index_column(self, index_name, table, column):
        """Called by interfaces to index specific column in table"""
        cursor = self.get_cursor()
        try:
            cursor.execute(
                'CREATE INDEX {0} on {1}({2})'.format(index_name, table, column))
        except sqlite3.OperationalError as error:
            print(error)
            print('Skipping index creation and assuming it exists already')
        else:
            self.conn.commit()

    def get_inclause(self, inlist):
        """Returns SQL IN clauses"""
        return 'IN ({0})'.format(', '.join('?' * len(inlist)))

    def get_sql_select(self, columns, table, distinct=False):
        """Creates and returns an SQL SELECT statement"""
        sql = 'SELECT {0} {1} FROM {2}'
        dist = {True: 'DISTINCT', False: ''}[distinct]
        return sql.format(dist, ', '.join(columns), table)

    def store_many(self, sql, values):
        """Abstraction over executemany method"""
        cursor = self.get_cursor()
        cursor.executemany(sql, values)
        self.conn.commit()

    def execute_sql(self, sql):
        """Executes SQL and returns cursor for it"""
        cursor = self.get_cursor()
        cursor.execute(sql)
        return cursor


class ResultLookupInterface(DatabaseConnection):
    """Connection subclass shared by result lookup interfaces that need
    to query the database when storing to get e.g. spectra id"""

    def get_mzmlfile_map(self):
        """Returns dict of mzmlfilenames and their db ids"""
        cursor = self.get_cursor()
        cursor.execute('SELECT mzmlfile_id, mzmlfilename FROM mzmlfiles')
        return {fn: fnid for fnid, fn in cursor.fetchall()}

    def get_spectra_id(self, fn_id, retention_time=None, scan_nr=None):
        """Returns spectra id for spectra filename and retention time"""
        cursor = self.get_cursor()
        sql = 'SELECT spectra_id FROM mzml WHERE mzmlfile_id=? '
        values = [fn_id]
        if retention_time is not None:
            sql = '{0} AND retention_time=?'.format(sql)
            values.append(retention_time)
        if scan_nr is not None:
            sql = '{0} AND scan_nr=?'.format(sql)
            values.append(scan_nr)
        cursor.execute(sql, tuple(values))
        return cursor.fetchone()[0]
