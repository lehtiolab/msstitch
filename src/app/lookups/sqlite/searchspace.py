from app.lookups.sqlite.base import DatabaseConnection


class SearchSpaceDB(DatabaseConnection):
    def add_tables(self, tabletypes):
        """Creates a searchspace lookup sqlite."""
        self.create_tables(['known_searchspace', 'protein_peptides', 'proteins', 'protein_seq'])

    def write_peps(self, peps):
        """Writes peps to db. We can reverse to be able to look up
        peptides that have some amino acids missing at the N-terminal.
        This way we can still use the index.
        """
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT OR IGNORE INTO known_searchspace(seqs) VALUES (?)', peps)
        self.conn.commit()

    def index_peps(self, reverse_seqs):
        if reverse_seqs:
            self.index_column('reverse_seqs_index', 'known_searchspace',
                              'seqs COLLATE NOCASE', unique=True)
        else:
            self.index_column('sequence_index', 'known_searchspace',
                              'seqs', unique=True)

    def store_pep_proteins(self, pepproteins):
        cursor = self.get_cursor()
        cursor.executemany('INSERT INTO proteins(protein_acc) VALUES(?)', ((x,) for x in pepproteins))
        cursor.executemany('INSERT INTO protein_seq(protein_acc, sequence) VALUES(?,?)',
                ((k, v['seq']) for k, v in pepproteins.items()))
        cursor.executemany('INSERT INTO protein_peptides(seq, protein_acc, pos) '
                           'VALUES(?, ?, ?)', ((pep[0], k, pep[1]) for k,v in pepproteins.items()
                               for pep in v['peps']))
        self.conn.commit()

    def store_tryp_peps_mapped(self, pepproteins):
        """Writes peps to db, but this time peps is a dict where {prot1: [PEPSEQ, IAMAPEPTIDE, ...], ...},
        and we also store the mapping to a peptide's proteins. We can reverse to be able to look up
        peptides that have some amino acids missing at the N-terminal.
        This way we can still use the index.
        """
        cursor = self.get_cursor()
        cursor.executemany('INSERT INTO proteins(protein_acc) VALUES(?)', ((x,) for x in pepproteins))
        cursor.executemany('INSERT INTO protein_peptides(seq, protein_acc, pos) VALUES(?, ?, ?)',
                ((pep, prot, 0) for prot,peps in pepproteins.items() for pep in peps))

    def index_proteins(self):
        self.index_column('pepix', 'protein_peptides', 'seq')
        self.index_column('pp_protix', 'protein_peptides', 'protein_acc')
        self.index_column('pseq_protix', 'protein_seq', 'protein_acc')
        self.conn.commit()

    def index_peps_mapped(self, reverse_seqs):
        collate = ' COLLATE NOCASE' if reverse_seqs else ''
        self.index_column('reverse_seqs_index', 'protein_peptides', f'seq{collate}', unique=True)
        self.index_column('pp_protix', 'protein_peptides', 'protein_acc')
        self.conn.commit()

    def get_multi_seq(self, allseqs):
        cursor = self.get_cursor()
        maxparam = 999
        allseqs_found = set()
        for i in range(0, len(allseqs), maxparam):
            seqs = allseqs[i:i+maxparam]
            sql = 'SELECT seqs FROM known_searchspace WHERE seqs IN ({})'.format(
                    ', '.join('?' for _ in seqs))
            [allseqs_found.add(x[0]) for x in cursor.execute(sql, seqs)]
        return allseqs_found

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

    def get_proteins_from_peps(self, peptides, minpeplen):
        '''Retrieve peptide/proteinseq combinations'''
        cursor = self.get_cursor()
        sql = ('SELECT pp.protein_acc, pp.pos, ps.sequence FROM protein_peptides AS pp '
            'INNER JOIN protein_seq AS ps ON ps.protein_acc=pp.protein_acc WHERE pp.seq=?')
        return {seq: [(protid, pos, pseq) for protid, pos, pseq in
            cursor.execute(sql, (seq[:minpeplen],))] for seq in peptides}
