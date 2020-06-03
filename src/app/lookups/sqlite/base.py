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
                            'scan_sid TEXT',
                            'charge INTEGER',
                            'mz DOUBLE',
                            'retention_time DOUBLE',
                            'ion_injection_time DOUBLE',
                            'FOREIGN KEY(mzmlfile_id)'
                            'REFERENCES mzmlfiles'],
                   'ioninjtime': ['spectra_id TEXT',
                                  'ion_injection_time DOUBLE',
                                  'FOREIGN KEY(spectra_id)'
                                  'REFERENCES mzml',
                                  ],
                   'ionmob': ['spectra_id TEXT',
                                  'ion_mobility DOUBLE',
                                  'FOREIGN KEY(spectra_id)'
                                  'REFERENCES mzml',
                                  ],
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
                   'ms1_fwhm': ['feature_id INTEGER PRIMARY KEY',
                                'fwhm REAL',
                                'FOREIGN KEY(feature_id)'
                                'REFERENCES ms1_quant'],
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
                    'gn_id INTEGER',
                    'genetable_id INTEGER',
                    'quant REAL',
                    'FOREIGN KEY(gn_id) '
                    'REFERENCES associated_ids(gn_id) '
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
                                         'gn_id INTEGER',
                                         'channel_id INTEGER',
                                         'quantvalue REAL',
                                         'amount_psms INTEGER',
                                         'FOREIGN KEY(gn_id) '
                                         'REFERENCES associated_ids(gn_id) '
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
                   'gene_fdr': ['gene_id INTEGER',
                                'genetable_id INTEGER',
                                'fdr DOUBLE',
                                'FOREIGN KEY(gene_id) '
                                'REFERENCES genes(gene_id) '
                                'FOREIGN KEY(genetable_id) '
                                'REFERENCES '
                                'gene_tables(genetable_id)'
                                ],
                   'assoc_fdr': ['gn_id INTEGER',
                                 'genetable_id INTEGER',
                                 'fdr DOUBLE',
                                 'FOREIGN KEY(gn_id) '
                                 'REFERENCES associated_ids(gn_id) '
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
                                            'pacc_id INTEGER ',
                                            'FOREIGN KEY(pacc_id) '
                                            'REFERENCES '
                                            'proteins(pacc_id)'],
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
                             'gene_acc TEXT'],
                   'associated_ids': ['gn_id INTEGER PRIMARY KEY',
                                      'assoc_id TEXT'],
                   'ensg_proteins': ['gene_id INTEGER',
                           'pacc_id INTEGER',
                           'FOREIGN KEY(gene_id) REFERENCES genes(gene_id) '
                           'FOREIGN KEY(pacc_id) REFERENCES proteins(pacc_id)'],
                   'genename_proteins': ['gn_id INTEGER',
                           'pacc_id INTEGER',
                           'FOREIGN KEY(gn_id) REFERENCES associated_ids(gn_id) '
                           'FOREIGN KEY(pacc_id) REFERENCES proteins(pacc_id)'],
                           
                   'prot_desc': ['pacc_id INTEGER',
                                 'description TEXT',
                                 'FOREIGN KEY(pacc_id) '
                                 'REFERENCES proteins(pacc_id)'],
                   'known_searchspace': ['seqs TEXT UNIQUE'],
                   'protein_peptides': ['seq TEXT', 'protid TEXT',
                                        'pos INTEGER'],
                   }



"""
PGM:
master, gene, assoc, pgc, nrpgc, cov (all not pooled), put it in a dict lookup?

SELECT pgm.protein_acc, g.gene_acc, aid.assoc_id, cov.coverage, sub.pgc, sub.pgcnr from protein_group_master as pgm left outer join (select master_id, group_concat(protein_acc) as pgc, count(protein_acc) as pgcnr  from protein_group_content group by master_id) as sub on sub.master_id=pgm.master_id left outer join protein_coverage as cov on pgm.protein_acc=cov.protein_acc left outer join genes as g on pgm.protein_acc=g.protein_acc left outer join associated_ids as aid on aid.protein_acc=pgm.protein_acc;

nr_psm, nrpep nrunipep: (FIXME add FDR, quant
    -----
    SELECT master_id, bs.set_id, COUNT(DISTINCT ppg.psm_id), COUNT (DISTINCT ps.pep_id), COUNT(DISTINCT uni.pep_id), 
        FROM psm_protein_groups AS ppg 
        INNER JOIN psms ON ppg.psm_id=psms.psm_id 
        INNER JOIN peptide_sequences AS ps ON psms.pep_id=ps.pep_id
        INNER JOIN mzml ON mzml.spectra_id=psms.spectra_id
        INNER JOIN mzmlfiles AS mzf ON mzf.mzmlfile_id=mzml.mzmlfile_id
        INNER JOIN biosets AS bs ON bs.set_id=mzf.set_id 
        LEFT OUTER JOIN (
                SELECT ppg.pep_id AS pep_id from (
                    SELECT psms.pep_id AS pep_id, COUNT (DISTINCT ppg.master_id) AS nrpep 
                        FROM psm_protein_groups AS ppg INNER JOIN psms USING(psm_id)
                        GROUP BY psms.pep_id
                    ) AS ppg WHERE ppg.nrpep==1
                ) AS uni ON uni.pep_id=ps.pep_id
        GROUP BY master_id, bs.set_id


    select * from (select master_id, bs.set_id, count(ppg.psm_id) as nrpsm, COUNT(distinct ps.pep_id) AS nrpep from psm_protein_groups AS ppg inner join psms ON ppg.psm_id=psms.psm_id inner join peptide_sequences AS ps on psms.pep_id=ps.pep_id inner join mzml on mzml.spectra_id=psms.spectra_id inner join mzmlfiles as mzf on mzf.mzmlfile_id=mzml.mzmlfile_id inner join biosets as bs on bs.set_id=mzf.set_id group by master_id, bs.set_id) as sub where sub.nrpsm>2 ;


nr_unipep, fdr, ms1, iso


select pgm.protein_acc from protein_group_master as pgm inner join (select psm_id from group_concat(protein_psm) group_by protein_acc) as pp on pgm.protein_acc=pp.protein_acc;

"""
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
                print('WARNING: Table {} already exists in database, will '
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
        cur.execute('PRAGMA synchronous=OFF')

    def get_cursor(self):
        """Quickly get cursor, abstracting connection"""
        return self.conn.cursor()

    def close_connection(self):
        """Close connection to db, abstracts connection object"""
        self.conn.close()

    def index_column(self, index_name, table, column, unique=False):
        """Called by interfaces to index specific column in table"""
        cursor = self.get_cursor()
        unique = 'unique' if unique else ''
        try:
            cursor.execute(
                'CREATE {3} INDEX {0} on {1}({2})'.format(index_name, table, column, unique))
        except sqlite3.OperationalError as error:
            # Skipping index creation and assuming it exists already
            pass
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
