from app.lookups.sqlite.base import ResultLookupInterface


class PSMDB(ResultLookupInterface):
    def add_tables(self):
        self.create_tables(['psms', 'psmrows', 'peptide_sequences',
                            'proteins', 'protein_evidence', 'protein_seq',
                            'prot_desc'])

    def store_proteins(self, proteins, evidence_lvls=False, sequences=False):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO proteins(protein_acc) '
            'VALUES(?)', proteins)
        self.conn.commit()
        cursor = self.get_cursor()
        if evidence_lvls:
            cursor.executemany(
                'INSERT INTO protein_evidence(protein_acc, evidence_lvl) '
                'VALUES(?, ?)', evidence_lvls)
        if sequences:
            cursor.executemany(
                'INSERT INTO protein_seq(protein_acc, sequence) '
                'VALUES(?, ?)', sequences)
        self.conn.commit()
        self.index_column('proteins_index', 'proteins', 'protein_acc')

    def store_descriptions(self, descriptions):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO prot_desc(protein_acc, description) '
            'VALUES(?, ?)', descriptions)
        self.conn.commit()
        self.index_column('protdesc_index', 'prot_desc', 'protein_acc')

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
