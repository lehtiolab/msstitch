from app.lookups.sqlite.base import ResultLookupInterface


class PSMDB(ResultLookupInterface):
    def add_tables(self):
        self.create_tables(['psms', 'psmrows', 'peptide_sequences',
                            'proteins', 'protein_evidence', 'protein_seq',
                            'prot_desc', 'protein_psm', 'genes',
                            'associated_ids'])

    def store_fasta(self, prot, evids, seq, desc, ensg, symbols):
        self.store_proteins(prot, evidence_lvls=evids, sequences=seq)
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO prot_desc(protein_acc, description) '
            'VALUES(?, ?)', desc)
        self.conn.commit()
        self.index_column('protdesc_index', 'prot_desc', 'protein_acc')
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO genes(gene_acc, protein_acc) VALUES(?, ?)', ensg)
        self.conn.commit()
        self.index_column('gene_index', 'genes', 'protein_acc')
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO associated_ids(assoc_id, protein_acc) VALUES(?, ?)',
            symbols)
        self.conn.commit()
        self.index_column('associd_index', 'associated_ids', 'protein_acc')

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
        self.index_column('evidence_index', 'protein_evidence', 'protein_acc')

    def get_protein_gene_map(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT p.protein_acc, g.gene_acc, aid.assoc_id, d.description '
            'FROM proteins AS p '
            'LEFT OUTER JOIN genes AS g ON p.protein_acc=g.protein_acc '
            'LEFT OUTER JOIN associated_ids AS aid '
            'ON p.protein_acc=aid.protein_acc '
            'LEFT OUTER JOIN prot_desc AS d ON p.protein_acc=d.protein_acc'
        )
        gpmap = {p_acc: {'gene': gene, 'symbol': sym, 'desc': desc}
                 for p_acc, gene, sym, desc in cursor}
        return gpmap

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
                                    psm['spec_id'])
                                   for psm in psms))
        cursor.executemany(
            'INSERT INTO psmrows(psm_id, rownr) VALUES(?, ?)',
            ((psm['psm_id'], psm['rownr']) for psm in psms))
        self.conn.commit()

    def store_peptides_proteins(self, allpepprot, psmids_to_store):
        ppmap = {psm_id: allpepprot[psm_id] for psm_id in psmids_to_store}
        prot_psm_ids = ((prot_acc, psm_id)
                        for psm_id, prots in ppmap.items()
                        for prot_acc in prots)
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO protein_psm(protein_acc, psm_id)'
            ' VALUES (?, ?)', prot_psm_ids)
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
        self.index_column('pepid_index', 'peptide_sequences', 'pep_id')
        self.index_column('psmspepid_index', 'psms', 'pep_id')

    def index_protein_peptides(self):
        self.index_column('protein_index', 'protein_psm', 'protein_acc')
        self.index_column('protpsmid_index', 'protein_psm', 'psm_id')
