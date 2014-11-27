import sqlite3
import os
from tempfile import mkstemp


class DatabaseConnection(object):
    def __init__(self, fn=None, foreign_keys=False):
        self.fn = fn
        if self.fn is not None:
            self.connect(self.fn)

    def get_fn(self):
        return self.fn

    def create_db(self, workdir, tables, outfn=None, foreign_keys=False):
        """Creates a sqlite db file.
        tables is a dict with keys=table names, values=lists of cols.
        """
        if outfn is None:
            fd, outfn = mkstemp(prefix='msstitcher_tmp_', dir=workdir)
            os.close(fd)
        self.fn = outfn
        self.connect(outfn, foreign_keys)
        self.create_tables(tables)

    def create_tables(self, tables):
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
