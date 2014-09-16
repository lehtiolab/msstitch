import sqlite3
import os
from tempfile import mkstemp


class DatabaseConnection(object):
    def __init__(self, fn=None):
        self.fn = fn
        if self.fn is not None:
            self.connect(self.fn)

    def create_db(self, tables, outfn=None):
        """Creates a sqlite db file.
        tables is a dict with keys=table names, values=lists of cols.
        """
        if outfn is None:
            fd, outfn = mkstemp(prefix='msstitcher_tmp_')
            os.close(fd)
        self.fn = outfn
        self.connect(outfn)
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
    def create_searchspacedb(self, outfn):
        self.create_db({'known_searchspace': ['seqs TEXT']}, outfn)

    def write_peps(self, peps, reverse_seqs):
        """Writes peps to db. We can reverse to be able to look up
        peptides that have some amino acids missing at the N-terminal.
        This way we can still use the index.
        """
        if reverse_seqs:
            peps = [(x[0][::-1],) for x in peps]
        self.conn.executemany(
            'INSERT INTO known_searchspace(seqs) VALUES (?)', peps)
        self.conn.commit()

    def index_peps(self):
        self.conn.execute('CREATE INDEX seqs_index ON known_searchspace(seqs)')
        self.conn.commit()

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
        cur = self.conn.execute(sql, (seq, ))
        return cur.fetchone()[0] == 1


class QuantDB(DatabaseConnection):
    def create_quantdb(self):
        self.create_db({'quant': ['spectra_filename TEXT', 'scan_nr TEXT',
                                  'quantmap TEXT', 'intensity REAL']})

    def store_quants(self, quants):
        self.conn.executemany(
            'INSERT INTO quant(spectra_filename, scan_nr, quantmap, intensity)'
            ' VALUES (?, ?, ?, ?)',
            quants)
        self.conn.commit()

    def index_quants(self):
        self.conn.execute('CREATE INDEX spec_index ON quant(spectra_filename)')
        self.conn.execute('CREATE INDEX scan_index ON quant(scan_nr)')
        self.conn.commit()

    def lookup_quant(self, spectrafile, scannr):
        cur = self.conn.execute(
            'SELECT quantmap, intensity FROM quant WHERE '
            'spectra_filename=? AND scan_nr=?', (spectrafile, scannr))
        return cur.fetchall()

    def get_all_quantmaps(self):
        cur = self.conn.execute(
            'SELECT DISTINCT quantmap FROM quant')
        return cur.fetchall()

class PeptideFilterDB(DatabaseConnection):
    pass
