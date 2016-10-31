from app.lookups.sqlite.base import DatabaseConnection


class SearchSpaceDB(DatabaseConnection):
    def add_tables(self):
        """Creates a searchspace lookup sqlite."""
        self.create_tables(['known_searchspace', 'protein_peptides'])

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

    def index_peps(self, reverse_seqs):
        if reverse_seqs:
            self.index_column('reverse_seqs_index', 'known_searchspace',
                              'seqs COLLATE NOCASE')
        else:
            self.index_column('reverse_seqs_index', 'known_searchspace',
                              'seqs')

    def store_pep_proteins(self, pepproteins):
        cursor = self.get_cursor()
        cursor.executemany('INSERT INTO protein_peptides(seq, protid, pos) '
                           'VALUES(?, ?, ?)', pepproteins)
        self.conn.commit()

    def index_proteins(self):
        self.index_column('pepix', 'protein_peptides', 'seq')
        self.conn.commit()

    def check_seq_exists(self, seq, amount_ntermwildcards):
        """Look up sequence in sqlite DB. Returns True or False if it
        exists (or not). When looking up a reversed DB with
        ntermwildcards: we reverse the sequence of the pep and add
        a LIKE and %-suffix to the query.
        """
        cursor = self.get_cursor()
        if amount_ntermwildcards > 0:
            seq = seq[::-1]
            sqlseq = '{}{}'.format(seq, '%')
            # FIXME non-parametrized string binding because if ? binding is
            # used the INDEX is not used when looking up, apparently because
            # the query cant be optimized when using LIKE and binding.
            sql = ('select seqs from known_searchspace where seqs LIKE '
                   '"{}"'.format(sqlseq))
            for match in cursor.execute(sql):
                if match[0][:-amount_ntermwildcards] in seq:
                    return True
            return False
        else:
            sql = ('select exists(select seqs from known_searchspace '
                   'where seqs=? limit 1)')
            return cursor.execute(sql, (seq, )).fetchone()[0] == 1

    def get_protein_from_pep(self, peptide):
        cursor = self.get_cursor()
        cursor.execute('SELECT protid, pos FROM protein_peptides WHERE seq='
                       '"{}"'.format(peptide))
        return cursor
