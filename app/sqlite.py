import sqlite3
import os
from tempfile import mkstemp


class DatabaseConnection(object):
    def __init__(self, fn=None, foreign_keys=False):
        self.fn = fn
        if self.fn is not None:
            self.connect(self.fn)

    def create_db(self, tables, outfn=None, foreign_keys=False):
        """Creates a sqlite db file.
        tables is a dict with keys=table names, values=lists of cols.
        """
        if outfn is None:
            fd, outfn = mkstemp(prefix='msstitcher_tmp_', dir=os.getcwd())
            os.close(fd)
        self.fn = outfn
        self.connect(outfn, foreign_keys)
        cursor = self.get_cursor()
        for table in tables:
            columns = tables[table]
            cursor.execute('CREATE TABLE {0}({1})'.format(
                table, ', '.join(columns)))
        self.conn.commit()

    def connect(self, fn, foreign_keys=False):
        self.conn = sqlite3.connect(fn)
        cur = self.get_cursor()
        if foreign_keys:
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
        cursor = self.get_cursor()
        cursor.executemany(
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
        cursor = self.get_cursor()
        cursor.execute(sql, (seq, ))
        return cursor.fetchone()[0] == 1


class QuantDB(DatabaseConnection):
    def create_quantdb(self):
        self.create_db({'quant': ['spectra_filename TEXT', 'scan_nr TEXT',
                                  'quantmap TEXT', 'intensity REAL']})

    def store_quants(self, quants):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO quant(spectra_filename, scan_nr, quantmap, intensity)'
            ' VALUES (?, ?, ?, ?)',
            quants)
        self.conn.commit()

    def index_quants(self):
        self.index_column('spec_index', 'quant', 'spectra_filename')
        self.index_column('scan_index', 'quant', 'scan_nr')

    def lookup_quant(self, spectrafile, scannr):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT quantmap, intensity FROM quant WHERE '
            'spectra_filename=? AND scan_nr=?', (spectrafile, scannr))
        return cursor.fetchall()

    def get_all_quantmaps(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT quantmap FROM quant')
        return cursor.fetchall()


MASTER_INDEX = 1
PROTEIN_ACC_INDEX = 2
PEPTIDE_COUNT_INDEX = 3
PSM_COUNT_INDEX = 4
PROTEIN_SCORE_INDEX = 5
EVIDENCE_LVL_INDEX = 6
COVERAGE_INDEX = 7


class ProteinGroupDB(DatabaseConnection):
    def create_pgdb(self):
        self.create_db({'psms': ['psm_id TEXT PRIMARY KEY NOT NULL',
                                 'sequence TEXT',
                                 'score TEXT'],
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
                                               ' proteins(protein_acc)']
                        }, foreign_keys=True)

    def store_proteins(self, proteins, evidence_lvls=False, sequences=False):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO proteins(protein_acc) '
            'VALUES(?)', proteins)
        if evidence_lvls:
            cursor.executemany(
                'INSERT INTO protein_evidence(protein_acc, evidence) '
                'VALUES(?, ?)', evidence_lvls)
        if sequences:
            cursor.executemany(
                'INSERT INTO protein_seq(protein_acc, sequence) '
                'VALUES(?, ?)', sequences)
        self.conn.commit()

    def store_peptides_proteins(self, ppmap, proteins_already_stored=False):
        def generate_proteins(pepprots):
            proteins = {}
            for psmvals in pepprots.values():
                for protein in psmvals['proteins']:
                    try:
                        proteins[protein]
                    except KeyError:
                        proteins[protein] = 1
                        yield (protein,)

        def generate_protein_psm_ids(pepprots):
            for psmvals in pepprots.values():
                for protein in psmvals['proteins']:
                    yield (protein, psmvals['psm_id'])

        def generate_psms(pepprots):
            for row, psmvals in pepprots.items():
                yield (psmvals['psm_id'], row, psmvals['seq'], psmvals['score'])

        psms = generate_psms(ppmap) 
        psms = sorted(psms, key=lambda x: x[1])  # sorts on psm rows
        prot_psm_ids = generate_protein_psm_ids(ppmap)
        if not proteins_already_stored:
            self.store_proteins(generate_proteins(ppmap))
        self.store_psm_relations(psms, prot_psm_ids)

    def store_psm_relations(self, psms, prot_psm_ids):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO psms(psm_id, sequence, score)'
            ' VALUES(?, ?, ?)', ((x[0], x[2], x[3]) for x in psms))
        cursor.executemany(
            'INSERT INTO psmrows(psm_id, rownr) VALUES(?, ?)',
            ((psm[0], psm[1]) for psm in psms))
        cursor.executemany(
            'INSERT INTO protein_psm(protein_acc, psm_id)'
            ' VALUES (?, ?)', prot_psm_ids)
        self.conn.commit()

    def index_protein_peptides(self):
        self.index_column('protein_index', 'protein_psm', 'protein_acc')
        self.index_column('psmid_index', 'protein_psm', 'psm_id')

    def store_masters(self, allmasters, psm_masters):
        allmasters = ((x,) for x in allmasters)
        psms = ((psm_id, master) for psm_id, masters in psm_masters.items()
                for master in masters)
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO protein_group_master(protein_acc) VALUES(?)',
            allmasters)
        cursor.executemany(
            'INSERT INTO psm_protein_groups(psm_id, master) '
            'VALUES(?, ?)', psms)
        self.conn.commit()

    def store_coverage(self, coverage):
        cursor = self.get_cursor()
        sql = ('INSERT INTO protein_coverage(protein_acc, coverage) '
               'VALUES(?, ?)')
        cursor.executemany(sql, coverage)
        self.conn.commit()

    def store_protein_group_content(self, protein_groups):
        cursor = self.get_cursor()
        cursor.executemany('INSERT INTO protein_group_content('
                           'protein_acc, master, peptide_count, '
                           'psm_count, protein_score) '
                           'VALUES(?, ?, ?, ?, ?)', protein_groups)
        self.conn.commit()

    def get_allpsms_masters(self):
        sql = ('SELECT pgm.protein_acc, pp.psm_id FROM protein_group_master AS pgm '
               'JOIN protein_psm as pp USING(protein_acc)')
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

    def get_peptides_from_protein(self, protein):
        pepsql = self.get_sql_select(['psm_id'], 'protein_psm',
                                     distinct=True)
        pepsql = '{0} WHERE protein_acc=?'.format(
            pepsql)
        cursor = self.get_cursor()
        peptides = cursor.execute(pepsql, (protein,)).fetchall()
        return [x[0] for x in peptides]

    def get_proteins_from_psms(self, psms):
        protsql = self.get_sql_select(['protein_acc'],
                                      'protein_psm', distinct=True)
        protsql = '{0} WHERE psm_id {1}'.format(
            protsql, self.get_inclause(psms))
        cursor = self.get_cursor()
        return [x[0] for x in cursor.execute(protsql, psms).fetchall()]

    def get_all_proteins_psms_seq(self):
        sql = ('SELECT p.protein_acc, ps.sequence, pp.psm_id, psms.sequence '
               'FROM proteins AS p '
               'JOIN protein_seq AS ps USING(protein_acc) '
               'JOIN protein_psm AS pp USING(protein_acc) '
               'JOIN psms AS psms USING(psm_id)'
               )
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_all_psms_proteingroups(self, coverage, evidence_levels):
        fields = ['pr.rownr', 'ppg.master', 'pgc.protein_acc',
                  'pgc.peptide_count', 'pgc.psm_count', 'pgc.protein_score']
        joins = [('psm_protein_groups', 'ppg', 'psm_id'),
                 ('protein_group_content', 'pgc', 'master')]
        if evidence_levels:
            fields.append('pev.evidence_lvl')
            joins.append(('protein_evidence', 'pev', 'protein_acc'))
        if coverage:
            fields.append('psq.sequence')
            joins.append(('protein_seq', 'psq', 'protein_acc'))
        join_sql = '\n'.join(['JOIN {0} AS {1} USING({2})'.format(
            j[0], j[1], j[2]) for j in joins])
        sql = 'SELECT {0} FROM psmrows AS pr {1}'.format(
            ', '.join(fields), join_sql)
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_proteins_peptides_from_psms(self, psms):
        """Returns dict of proteins and lists of corresponding peptides
        from db. DB call gets all rows where psm_id is in peptides.
        """
        sql = ('SELECT protein_psm.protein_acc, psms.sequence, psms.score, '
               'psm_id FROM psms '
               'JOIN protein_psm USING(psm_id)')
        sql = '{0} WHERE psm_id {1}'.format(sql, self.get_inclause(psms))
        cursor = self.get_cursor()
        return cursor.execute(sql, psms).fetchall()

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
        cursor = self.get_cursor()
        proteins_not_in_group = cursor.execute(not_in_sql,
                                               proteins + peptides)
        return [x[0] for x in proteins_not_in_group]
