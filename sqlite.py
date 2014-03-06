import sqlite3
import os
from tempfile import mkstemp


class DatabaseConnection(object):
    def __init__(self, fn=None):
        self.fn = fn

    def create_db(self, tables):
        """Creates a sqlite db file.
        tables is a dict with keys=table names, values=lists of cols.
        """
        fd, self.fn = mkstemp(prefix='pycolator_tmp_')
        os.close(fd)
        self.connect_searchspace(self.fn)
        for table in tables:
            columns = tables[table]
            sql = 'CREATE TABLE ?({0})'.format(','.join(['?'] * len(columns)))
            sql_names = (table, ) + tuple(columns)
            self.conn.execute(sql, sql_names)
        self.conn.commit()

    def connect_searchspace(self, fn):
        self.conn = sqlite3.connect(fn)

    def close_connection(self):
        self.conn.close()

    def create_searchspacedb(self):
        self.create_db({'known_searchspace': ['seqs']})

    def write_peps(self, peps):
        self.conn.executemany(
            'INSERT INTO known_searchspace(seqs) VALUES (?)', peps)
        self.conn.commit()

    def index_peps(self):
        self.conn.execute('CREATE INDEX seqs_index ON known_searchspace(seqs)')
        self.conn.commit()

    def check_seq_exists(self, seq):
        cur = self.conn.execute(
            'select exists(select seqs from known_searchspace '
            'where seqs=? limit 1)',
            (seq, ))
        return cur.fetchone()[0] == 1
