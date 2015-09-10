from app.lookups.sqlite.base import ResultLookupInterface


class PSMDB(ResultLookupInterface):
    def add_tables(self):
        self.create_tables(['psms', 'psmrows', 'peptide_sequences'])

    def store_pepseqs(self, sequences):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO peptide_sequences(sequence) VALUES(?)', sequences)
        self.conn.commit()

    def store_psms(self, psms):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO psms(psm_id, pep_id, score, spectra_id) '
            'VALUES(?, ?, ?, ?)', ((psm['psm_id'], psm['seq'], psm['score'],
                                    self.get_spectra_id(psm['specfn'],
                                                        scan_nr=psm['scannr']))
                                   for psm in psms))
        cursor.executemany(
            'INSERT INTO psmrows(psm_id, rownr) VALUES(?, ?)',
            ((psm['psm_id'], psm['rownr']) for psm in psms))
        self.conn.commit()

    def get_peptide_seq_map(self):
        cursor = self.get_cursor()
        seqs = cursor.execute('SELECT pep_id, sequence FROM peptide_sequences')
        return {seq: pepid for pepid, seq in seqs.fetchall()}

    def index_psms(self):
        self.index_column('psmid_index', 'psms', 'psm_id')
        self.index_column('psmspecid_index', 'psms', 'spectra_id')
        self.index_column('psmrowid_index', 'psmrows', 'psm_id')
        self.index_column('psmrow_index', 'psmrows', 'rownr')
        self.index_column('pepseq_index', 'peptide_sequences', 'sequence')
