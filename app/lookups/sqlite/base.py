import sqlite3
from app.lookups.sqlite import proteingroups, quant, searchspace, \
    biosets, spectra


mslookup_tables = {'biosets': ['set_id INTEGER PRIMARY KEY',
                               'set_name TEXT'],
                   'mzmlfiles': ['mzmlfile_id INTEGER PRIMARY KEY',
                                 'mzmlfilename TEXT',
                                 'set_id INTEGER',
                                 'FOREIGN KEY(set_id)'
                                 'REFERENCES biosets'
                                 ],
                   'mzml': ['spectra_id INTEGER PRIMARY KEY',
                            'mzmlfile_id INTEGER',
                            'scan_nr TEXT',
                            'retention_time REAL',
                            'FOREIGN KEY(mzmlfile_id)'
                            'REFERENCES mzmlfiles'],
                   'isobaric_quant': ['mzmlfile_id INTEGER',
                                      'retention_time REAL',
                                      'quantmap TEXT',
                                      'intensity REAL',
                                      'FOREIGN KEY(mzmlfile_id)'
                                      'REFERENCES mzmlfiles '
                                      'FOREIGN KEY(retention_time)'
                                      'REFERENCES mzml'
                                      ],
                   'ms1_quant': ['mzmlfile_id INTEGER',
                                 'retention_time REAL', 'mz REAL',
                                 'charge INTEGER', 'intensity REAL',
                                 'FOREIGN KEY(mzmlfile_id)'
                                 'REFERENCES mzmlfiles '
                                 'FOREIGN KEY(retention_time)'
                                 'REFERENCES mzml'],
                   'psms': ['psm_id TEXT PRIMARY KEY NOT NULL',
                            'sequence TEXT',
                            'score TEXT',
                            'set_name TEXT',
                            'spectra_id INTEGER',
                            'FOREIGN KEY(spectra_id)'
                            'REFERENCES mzml'
                            ],
                   # FIXME create a spectra_id: scan nr reference when
                   # inserting psms, otherwise we must look them up?
                   # protein table - protein ref, quant data, set ref
                   'psmrows': ['psm_id TEXT',
                               'rownr INTEGER',
                               'FOREIGN KEY(psm_id) '
                               'REFERENCES psms(psm_id)'],
                   'proteins': ['protein_acc TEXT PRIMARY KEY NOT NULL'],
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
                   'protein_group_master': ['protein_acc TEXT',
                                            'FOREIGN KEY(protein_acc) '
                                            'REFERENCES '
                                            'proteins(protein_acc)'],
                   'protein_group_content': ['protein_acc TEXT',
                                             'master TEXT',
                                             'peptide_count INTEGER',
                                             'psm_count INTEGER',
                                             'protein_score INTEGER',
                                             'FOREIGN KEY(protein_acc) '
                                             'REFERENCES '
                                             'proteins(protein_acc) '
                                             'FOREIGN KEY(master) '
                                             'REFERENCES '
                                             'proteins'
                                             '(protein_acc)'
                                             ],
                   'psm_protein_groups': ['psm_id TEXT',
                                          'master TEXT',
                                          'FOREIGN KEY(psm_id) REFERENCES'
                                          ' psms(psm_id)',
                                          'FOREIGN KEY(master) REFERENCES'
                                          ' proteins(protein_acc)'],
                   'prot_desc': ['protein_acc TEXT',
                                 'description TEXT',
                                 'FOREIGN KEY(protein_acc) '
                                 'REFERENCES proteins(protein_acc)'],
                   }


class DatabaseConnection(object):
    def __init__(self, fn=None):
        self.fn = fn
        if self.fn is not None:
            self.connect(self.fn)

    def get_fn(self):
        return self.fn

    def create_tables(self, tables):
        cursor = self.get_cursor()
        for table in tables:
            columns = mslookup_tables[table]
            cursor.execute('CREATE TABLE {0}({1})'.format(
                table, ', '.join(columns)))
        self.conn.commit()

    def connect(self, fn):
        self.conn = sqlite3.connect(fn)
        cur = self.get_cursor()
        cur.execute('PRAGMA FOREIGN_KEYS=ON')

    def get_cursor(self):
        return self.conn.cursor()

    def close_connection(self):
        self.conn.close()

    def index_column(self, index_name, table, column):
        cursor = self.get_cursor()
        cursor.execute(
            'CREATE INDEX {0} on {1}({2})'.format(index_name, table, column))
        self.conn.commit()

    def get_inclause(self, inlist):
        return 'IN ({0})'.format(', '.join('?' * len(inlist)))

    def get_sql_select(self, columns, table, distinct=False):
        sql = 'SELECT {0} {1} FROM {2}'
        dist = {True: 'DISTINCT', False: ''}[distinct]
        return sql.format(dist, ', '.join(columns), table)

    def store_many(self, sql, values):
        cursor = self.get_cursor()
        cursor.executemany(sql, values)
        self.conn.commit()


class ResultLookupInterface(DatabaseConnection):
    def get_mzmlfile_map(self):
        cursor = self.get_cursor()
        cursor.execute('SELECT mzmlfile_id, mzmlfilename FROM mzmlfiles')
        return cursor.fetchall()


def get_lookup(fn, lookuptype):
    lookupmap = {'biosets': biosets.BioSetDB,
                 'spectra': spectra.SpectraDB,
                 'proteingroups': proteingroups.ProteinGroupDB,
                 'quant': quant.QuantDB,
                 'searchspace': searchspace.SearchSpaceDB,
                 }
    return lookupmap[lookuptype](fn)


def create_new_lookup(fn, lookuptype):
    with open(fn, 'w'):
        pass
    return get_lookup(fn, lookuptype)
