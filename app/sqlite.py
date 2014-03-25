import sqlite3
import os
from tempfile import mkstemp


class DatabaseConnection(object):
    def __init__(self, fn=None):
        self.fn = fn
        if self.fn is not None:
            self.connect(self.fn)

    def create_db(self, tables):
        """Creates a sqlite db file.
        tables is a dict with keys=table names, values=lists of cols.
        """
        fd, self.fn = mkstemp(prefix='msstitcher_tmp_')
        os.close(fd)
        self.connect(self.fn)
        for table in tables:
            columns = tables[table]
            self.conn.execute('CREATE TABLE {0}({1})'.format(
                table, ', '.join(columns)))
        self.conn.commit()

    def connect(self, fn):
        self.conn = sqlite3.connect(fn)

    def close_connection(self):
        self.conn.close()


class SearchSpaceDB(DatabaseConnection):
    def create_searchspacedb(self):
        self.create_db({'known_searchspace': ['seqs TEXT']})

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


class QuantDB(DatabaseConnection):
    def create_quantdb(self):
        self.create_db({'quant': ['spectra_filename TEXT', 'scan_nr TEXT',
                                  'quantmap TEXT', 'intensity REAL']})

    def store_quants(self, quants):
        self.conn.executemany(
            'INSERT INTO quant(spectra_filename, scan_nr, quantmap, intensity)'
            ' VALUES (?, ?, ?, ?)',
            (quants,))

    def lookup_quant(self, spectrafile, scannr):
        cur = self.conn.execute(
            'SELECT quantmap, intensity FROM quant WHERE '
            'spectra_filename=? AND scan_nr=?', (spectrafile, scannr))
        return cur.fetchall()


class PeptideFilterDB(DatabaseConnection):
    pass
