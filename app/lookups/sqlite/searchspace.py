from app.lookups.sqlite.base import DatabaseConnection


class SearchSpaceDB(DatabaseConnection):
    def create_searchspacedb(self, outfn):
        """Creates a searchspace lookup sqlite. Since the ultimate output
        of this is the sqlite file, we use None as workdir."""
        self.create_db(None,
                       {'known_searchspace': ['seqs TEXT']}, outfn)

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
