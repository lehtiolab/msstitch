import sqlite3
import os
from tempfile import mkstemp


class DatabaseConnection(object):
    def __init__(self, fn=None, foreign_keys=False):
        self.fn = fn
        if self.fn is not None:
            self.cursor = self.connect(self.fn)

    def create_db(self, tables, outfn=None, foreign_keys=False):
        """Creates a sqlite db file.
        tables is a dict with keys=table names, values=lists of cols.
        """
        if outfn is None:
            fd, outfn = mkstemp(prefix='msstitcher_tmp_', dir=os.getcwd())
            os.close(fd)
        self.fn = outfn
        self.cursor = self.connect(outfn, foreign_keys)
        for table in tables:
            columns = tables[table]
            self.cursor.execute('CREATE TABLE {0}({1})'.format(
                table, ', '.join(columns)))
        self.conn.commit()

    def connect(self, fn, foreign_keys=False):
        self.conn = sqlite3.connect(fn)
        cur = self.conn.cursor()
        if foreign_keys:
            cur.execute('PRAGMA FOREIGN_KEYS=ON')
        return cur

    def close_connection(self):
        self.conn.close()

    def index_column(self, index_name, table, column):
        self.cursor.execute(
            'CREATE INDEX {0} on {1}({2})'.format(index_name, table, column))
        self.conn.commit()

    def get_inclause(self, inlist):
        return 'IN ({0})'.format(', '.join('?' * len(inlist)))

    def get_sql_select(self, columns, table, distinct=False):
        sql = 'SELECT {0} {1} FROM {2}'
        dist = {True: 'DISTINCT', False: ''}[distinct]
        return sql.format(dist, ', '.join(columns), table)


class SearchSpaceDB(DatabaseConnection):
    def create_searchspacedb(self, outfn):
        self.create_db({'known_searchspace': ['seqs TEXT']}, outfn)

    def write_peps(self, peps, reverse_seqs):
        """Writes peps to db. We can reverse to be able to look up
        peptides that have some amino acids missing at the N-terminal.
        This way we can still use the index.
        """
        if reverse_seqs:
            peps = [(x[0][::-1],) for x in peps]
        self.cursor.executemany(
            'INSERT INTO known_searchspace(seqs) VALUES (?)', peps)
        self.conn.commit()

    def index_peps(self):
        self.index_column('seqs_index', 'known_searchspace', 'seqs')

    def check_seq_exists(self, seq, ntermwildcards=False):
        """Look up sequence in sqlite DB. Returns True or False if it
        exists (or not). When looking up a reversed DB with
        ntermwildcards: we reverse the sequence of the pep and add
        a LIKE and %-suffix to the query.
        """
        if ntermwildcards:
            seq = seq[::-1]
            comparator, seqmod = ' LIKE ', '%'
        else:
            comparator, seqmod = '=', ''

        sql = ('select exists(select seqs from known_searchspace '
               'where seqs{0}? limit 1)'.format(comparator))
        seq = '{0}{1}'.format(seq, seqmod)
        self.cursor.execute(sql, (seq, ))
        return self.cursor.fetchone()[0] == 1


class QuantDB(DatabaseConnection):
    def create_quantdb(self):
        self.create_db({'quant': ['spectra_filename TEXT', 'scan_nr TEXT',
                                  'quantmap TEXT', 'intensity REAL']})

    def store_quants(self, quants):
        self.cursor.executemany(
            'INSERT INTO quant(spectra_filename, scan_nr, quantmap, intensity)'
            ' VALUES (?, ?, ?, ?)',
            quants)
        self.conn.commit()

    def index_quants(self):
        self.index_column('spec_index', 'quant', 'spectra_filename')
        self.index_column('scan_index', 'quant', 'scan_nr')

    def lookup_quant(self, spectrafile, scannr):
        self.cursor.execute(
            'SELECT quantmap, intensity FROM quant WHERE '
            'spectra_filename=? AND scan_nr=?', (spectrafile, scannr))
        return self.cursor.fetchall()

    def get_all_quantmaps(self):
        self.cursor.execute(
            'SELECT DISTINCT quantmap FROM quant')
        return self.cursor.fetchall()


MASTER_INDEX = 1
PROTEIN_ACC_INDEX = 2
PEPTIDE_COUNT_INDEX = 3
PSM_COUNT_INDEX = 4
PROTEIN_SCORE_INDEX = 5
EVIDENCE_LVL_INDEX = 6
COVERAGE_INDEX = 7


class ProteinGroupDB(DatabaseConnection):
    def create_pgdb(self):
        self.create_db({'psms': ['psm_id INTEGER PRIMARY KEY NOT NULL',
                                 'sequence TEXT',
                                 'score TEXT'],
                        'proteins': ['protein_acc TEXT'],
                        'protein_psm': ['protein_acc TEXT',
                                        'psm_id INTEGER',
                                        'FOREIGN KEY(protein_acc) '
                                        'REFERENCES proteins(protein_acc)',
                                        'FOREIGN KEY(psm_id) '
                                        'REFERENCES psms(psm_id)'],
                        'protein_evidence': ['protein_acc TEXT ',
                                             'evidence_lvl REAL '
                                             'FOREIGN KEY(protein_acc) '
                                             'REFERENCES '
                                             'proteins(protein_acc)'],
                        'protein_seq': ['protein_acc TEXT ',
                                        'sequence TEXT '
                                        'FOREIGN KEY(protein_acc) '
                                        'REFERENCES '
                                        'proteins(protein_acc)'],
                        'protein_group_master': ['master '
                                                 'TEXT PRIMARY KEY NOT NULL'],
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
                                                  'protein_group_master'
                                                  '(master)'
                                                  ],
                        'psm_protein_groups': ['psm_id INTEGER',
                                               'master TEXT',
                                               'FOREIGN KEY(psm_id) REFERENCES'
                                               ' psms(psm_id)',
                                               'FOREIGN KEY(master) REFERENCES'
                                               ' protein_group_master(master)']
                        }, foreign_keys=True)

    def store_peptides_proteins(self, ppmap):
        def generate_proteins(pepprots, only_proteins=False):
            for psm_id, psmvals in pepprots.items():
                for protein in psmvals['proteins']:
                    protein_out = {False: (protein, psm_id),
                                   True: (protein,)}[only_proteins]
                    yield protein_out
        psms = ((k, v['seq'], v['score']) for k, v in ppmap.items())
        self.cursor.executemany(
            'INSERT INTO psms(psm_id, sequence, score)'
            ' VALUES(?, ?, ?)', psms)
        self.cursor.executemany(
            'INSERT INTO proteins(protein_acc) VALUES(?)',
            generate_proteins(ppmap, True))
        self.cursor.executemany(
            'INSERT INTO protein_psm(protein_acc, psm_id)'
            ' VALUES (?, ?)', generate_proteins(ppmap))
        self.conn.commit()

    def index_protein_peptides(self):
        self.index_column('protein_index', 'protein_psm', 'protein_acc')

    def store_masters(self, allmasters, psms):
        print('Storing {0} masters for {1} PSMs'.format(len(allmasters), len(psms)))
        allmasters = ((x,) for x in allmasters)
        self.cursor.executemany(
            'INSERT INTO protein_group_master(master) VALUES(?)',
            allmasters)
        self.cursor.executemany(
            'INSERT INTO psm_protein_groups(psm_id, master) '
            'VALUES(?, ?)', psms)
        self.conn.commit()

    def store_protein_group_content(self, protein_groups):
        self.cursor.executemany('INSERT INTO protein_group_content('
                                'protein_acc, master, peptide_count, '
                                'psm_count, protein_score) '
                                'VALUES(?, ?, ?, ?, ?)', protein_groups)
        self.conn.commit()

    def get_all_masters(self):
        sql = self.get_sql_select(['master'], 'protein_group_master')
        return self.cursor.execute(sql).fetchall()

    def get_proteins_for_peptide(self, psm_id):
        """Returns list of proteins for a passed psm_id"""
        protsql = self.get_sql_select(['protein_acc'], 'protein_psm')
        protsql = '{0} WHERE psm_id=?'.format(protsql)
        proteins = self.cursor.execute(protsql, psm_id).fetchall()
        return [x[0] for x in proteins]

    def get_protpepmap_from_proteins(self, proteins):
        pepsql = self.get_sql_select(['protein_acc', 'psm_id'],
                                     'protein_psm',
                                     distinct=True)
        pepsql = '{0} WHERE protein_acc {1}'.format(
            pepsql, self.get_inclause(proteins))
        protpeps = self.cursor.execute(pepsql, proteins).fetchall()
        outmap = {}
        for protein, peptide in protpeps:
            try:
                outmap[protein].append(peptide)
            except KeyError:
                outmap[protein] = [peptide]
        return outmap

    def get_peptides_from_protein(self, protein):
        pepsql = self.get_sql_select(['psm_id'], 'protein_psm',
                                     distinct=True)
        pepsql = '{0} WHERE protein_acc=?'.format(
            pepsql)
        peptides = self.cursor.execute(pepsql, (protein,)).fetchall()
        return [x[0] for x in peptides]

    def get_proteins_from_psms(self, psms):
        protsql = self.get_sql_select(['protein_acc'],
                                      'protein_psm', distinct=True)
        protsql = '{0} WHERE psm_id {1}'.format(
            protsql, self.get_inclause(psms))
        return [x[0] for x in self.cursor.execute(protsql, psms).fetchall()]

    def get_all_psms_proteingroups(self, evidence_levels, fasta):
        fields = ['p.psm_id', 'ppg.master', 'pgc.protein_acc',
                  'pgc.peptide_count', 'pgc.psm_count', 'pgc.protein_score']
        joins = [('psms', 'p', 'psm_id'),
                 ('protein_group_content', 'pgc', 'master')]
        if evidence_levels:
            fields.append('pev.evidence_lvl')
            joins.append(('protein_evidence', 'pev', 'protein_acc'))
        if fasta:
            fields.append('psq.sequence')
            joins.append(('protein_seq', 'psq', 'protein_acc'))
        join_sql = '\n'.join(['JOIN {0} AS {1} USING({2})'.format(
            j[0], j[1], j[2]) for j in joins])
        sql = 'SELECT {0} FROM psm_protein_groups AS ppg {1}'.format(
            ', '.join(fields), join_sql)
        return self.cursor.execute(sql)

    def get_proteins_peptides_from_psms(self, psms):
        """Returns dict of proteins and lists of corresponding peptides
        from db. DB call gets all rows where psm_id is in peptides.
        """
        sql = ('SELECT protein_psm.protein_acc, psms.sequence, psms.score, '
               'psm_id FROM psms '
               'JOIN protein_psm USING(psm_id)')
        sql = '{0} WHERE psm_id {1}'.format(sql, self.get_inclause(psms))
        return self.cursor.execute(sql, psms).fetchall()

    def filter_proteins_with_missing_peptides(self, proteins, peptides):
        """Returns proteins of passed list that have peptides not in
        peptide list.
        """
        not_in_sql = self.get_sql_select(['protein_acc'], 'protein_psm',
                                         distinct=True)
        not_in_sql = '{0} WHERE protein_acc {1} AND psm_id NOT {2}'.format(
            not_in_sql,
            self.get_inclause(proteins),
            self.get_inclause(peptides))
        proteins_not_in_group = self.cursor.execute(not_in_sql,
                                                    proteins + peptides)
        return [x[0] for x in proteins_not_in_group]
