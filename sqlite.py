import sqlite3
import os
from tempfile import mkstemp


class DatabaseConnection(object):
    def __init__(self, fn=None):
        self.fn = fn
    
    def create_db(self):
        fd, self.fn = mkstemp(prefix='pycolator_tmp_', dir='/tmp')
        os.close(fd)
        self.conn = sqlite3.connect(self.fn)
        self.conn.execute("CREATE TABLE pycolator(seqs)")
        self.conn.commit()

    def connect_searchspace(self, fn):
        self.conn = sqlite3.connect(fn)

    def close_connection(self):
        self.conn.close()

    def write_peps(self, peps):
        self.conn.executemany("INSERT INTO pycolator(seqs) VALUES (?)", peps)
        self.conn.commit()

    def index_peps(self):
        self.conn.execute("CREATE INDEX seqs_index ON pycolator(seqs)")
        self.conn.commit()
    
    def check_seq_exists(self, seq):
        cur = self.conn.execute(
                "select exists(select seqs from pycolator where seqs=? limit 1)", 
                (seq,) )
        return cur.fetchone()[0] == 1

