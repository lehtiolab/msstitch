from app.lookups.sqlite.base import ResultLookupInterface


# Indices that belong to positions of these features in output from
# function get_all_psms_proteingroups:
MASTER_INDEX = 1
PROTEIN_ACC_INDEX = 2
PEPTIDE_COUNT_INDEX = 3
PSM_COUNT_INDEX = 4
PROTEIN_SCORE_INDEX = 5
COVERAGE_INDEX = 6
EVIDENCE_LVL_INDEX = 7


class PSMDB(ResultLookupInterface):
    def add_tables(self, tabletypes):
        self.create_tables(['psms', 'psmrows', 'peptide_sequences', 'fastafn',
                            'proteins', 'protein_evidence', 'protein_seq',
                            'prot_desc', 'protein_psm', 'genes',
                            'associated_ids', 'ensg_proteins', 'genename_proteins'])
        if 'proteingroup' in tabletypes:
            self.create_tables(['protein_coverage', 'protein_group_master',
                                'protein_group_content', 'psm_protein_groups'])

    def get_fasta_md5(self):
        cursor = self.get_cursor()
        cursor.execute('SELECT md5 FROM fastafn LIMIT 1')
        md5 = cursor.fetchone()
        return md5 if md5 is not None else False

    def store_fasta(self, fn, md5, prot, evids, seq, desc, ensg, symbols):
        cursor = self.get_cursor()
        cursor.execute('INSERT INTO fastafn(filename, md5) VALUES(?, ?)',
                (fn, md5))
        self.store_proteins(prot, evidence_lvls=evids, sequences=seq)
        protids = self.get_protids()
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO prot_desc(pacc_id, description) '
            'VALUES(?, ?)', ((protids[x[0]], x[1]) for x in desc))
        self.conn.commit()
        self.index_column('protdesc_index', 'prot_desc', 'pacc_id')

        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO genes(gene_acc) VALUES(?)', set((x[0],) for x in ensg))
        self.conn.commit()
        ensgids = {x[0]: x[1] for x in cursor.execute(
            'SELECT g.gene_acc, g.gene_id FROM genes AS g')}
        ensg = [(ensgids[x[0]], protids[x[1]]) for x in ensg]
        cursor.executemany(
            'INSERT INTO ensg_proteins(gene_id, pacc_id) VALUES(?, ?)', ensg)
        self.conn.commit()
        self.index_column('ensg_p_ix', 'ensg_proteins', 'pacc_id')
        self.index_column('ensg_ensg_ix', 'ensg_proteins', 'gene_id')

        cursor = self.get_cursor()
        # symbols is list of sym/prot id, so unique the symbol ids and save
        cursor.executemany(
            'INSERT INTO associated_ids(assoc_id) VALUES(?)',
            set((x[0],) for x in symbols))
        symids = {x[0]: x[1] for x in cursor.execute(
            'SELECT s.assoc_id, s.gn_id FROM associated_ids AS s')}
        symbols = [(symids[x[0]], protids[x[1]]) for x in symbols]
        cursor.executemany(
            'INSERT INTO genename_proteins(gn_id, pacc_id) VALUES(?, ?)',
            symbols)
        self.conn.commit()
        self.index_column('gp_p_ix', 'genename_proteins', 'pacc_id')
        self.index_column('gp_gn_ix', 'genename_proteins', 'gn_id')

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

    def get_protids(self):
        cur = self.get_cursor()
        return {x[0]: x[1] for x in cur.execute(
            'SELECT p.protein_acc, p.pacc_id FROM proteins AS p')}

    def get_protein_gene_map(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT p.protein_acc, g.gene_acc, aid.assoc_id, d.description '
            'FROM proteins AS p '
            'LEFT OUTER JOIN ensg_proteins AS ep ON p.pacc_id=ep.pacc_id '
            'LEFT OUTER JOIN genename_proteins AS gnp ON p.pacc_id=gnp.pacc_id '
            'LEFT OUTER JOIN genes AS g ON ep.gene_id=g.gene_id '
            'LEFT OUTER JOIN associated_ids AS aid ON gnp.gn_id=aid.gn_id '
            'LEFT OUTER JOIN prot_desc AS d ON p.pacc_id=d.pacc_id'
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
    
    def get_highest_rownr(self):
        cursor = self.get_cursor()
        cursor.execute('SELECT MAX(rownr) FROM psmrows')
        return int(cursor.fetchone()[0])

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
        self.index_column('proteinpsm_index', 'protein_psm', 'protein_acc')
        self.index_column('protpsmid_index', 'protein_psm', 'psm_id')

    def get_exp_spectra_data_rows(self, shiftrows):
        cursor = self.get_cursor()
        rowlim = 'WHERE pr.rownr>{} '.format(shiftrows - 1) if shiftrows else ''
        return cursor.execute('SELECT pr.rownr, bs.set_name, sp.retention_time, '
                              'iit.ion_injection_time, im.ion_mobility '
                              'FROM psmrows AS pr '
                              'JOIN psms AS p USING(psm_id) '
                              'JOIN mzml AS sp USING(spectra_id) '
                              'LEFT OUTER JOIN ioninjtime AS iit USING(spectra_id) '
                              'LEFT OUTER JOIN ionmob AS im USING(spectra_id) '
                              'JOIN mzmlfiles as mf USING(mzmlfile_id) '
                              'JOIN biosets AS bs USING(set_id) ' + rowlim +
                              'ORDER BY pr.rownr')

    def drop_psm_indices(self):
        cursor = self.get_cursor()
        for index in ['psmid_index', 'psmspecid_index', 'psmrowid_index',
                'psmrow_index', 'pepseq_index', 'pepid_index', 'psmspepid_index',
                'proteinpsm_index', 'protpsmid_index']:
            cursor.execute('DROP INDEX IF EXISTS {}'.format(index))
        self.conn.commit()

    def drop_pgroup_tables(self):
        for table in ['protein_coverage',
                'protein_group_content',
                'psm_protein_groups',
                'protein_group_master']:
            self.conn.execute('DROP TABLE IF EXISTS {}'.format(table))
        self.conn.commit()

    def store_masters(self, allmasters, psm_masters):
        protids = self.get_protids()
        allmasters = ((protids[x],) for x in allmasters)
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO protein_group_master(pacc_id) VALUES(?)',
            allmasters)
        master_ids = self.get_master_ids()
        psms = ((psm_id, master_ids[protids[master]])
                for psm_id, masters in psm_masters.items()
                for master in masters)
        cursor.executemany(
            'INSERT INTO psm_protein_groups(psm_id, master_id) '
            'VALUES(?, ?)', psms)
        self.conn.commit()
        self.index_column('psm_pg_index', 'psm_protein_groups', 'master_id')
        self.index_column('master_pacc_ix', 'protein_group_master', 'pacc_id')
        self.index_column('psm_pg_psmid_index', 'psm_protein_groups', 'psm_id')

    def update_master_proteins(self, new_masters):
        protids = self.get_protids()
        new_masters = ((protids[x[0]], x[1]) for x in new_masters)
        cur = self.get_cursor()
        sql = 'UPDATE protein_group_master SET pacc_id=? WHERE master_id=?'
        cur.executemany(sql, new_masters)
        self.conn.commit()

    def get_master_ids(self):
        cur = self.get_cursor()
        cur.execute('SELECT pacc_id, master_id FROM protein_group_master')
        return {p_acc: master_id for (p_acc, master_id) in cur}

    def store_coverage(self, coverage):
        cursor = self.get_cursor()
        sql = ('INSERT INTO protein_coverage(protein_acc, coverage) '
               'VALUES(?, ?)')
        cursor.executemany(sql, coverage)
        self.conn.commit()
        self.index_column('cov_index', 'protein_coverage', 'protein_acc')

    def store_protein_group_content(self, protein_groups):
        cursor = self.get_cursor()
        cursor.executemany('INSERT INTO protein_group_content('
                           'protein_acc, master_id, peptide_count, '
                           'psm_count, protein_score) '
                           'VALUES(?, ?, ?, ?, ?)', protein_groups)
        self.conn.commit()

    def index_protein_group_content(self):
        self.index_column('pgc_master_index', 'protein_group_content',
                          'master_id')

    def get_all_psm_protein_relations(self):
        sql = 'SELECT psm_id, protein_acc FROM protein_psm'
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_proteins_for_peptide(self, psm_id):
        """Returns list of proteins for a passed psm_id"""
        protsql = self.get_sql_select(['protein_acc'], 'protein_psm')
        protsql = '{0} WHERE psm_id=?'.format(protsql)
        cursor = self.get_cursor()
        proteins = cursor.execute(protsql, psm_id).fetchall()
        return [x[0] for x in proteins]

    def get_psms_for_proteins(self, proteins):
        sql = """SELECT protein_acc, psm_id FROM protein_psm 
        WHERE protein_acc IN ({})""".format(','.join('?' * len(proteins)))
        return self.conn.execute(sql, proteins)

    def get_all_proteins_psms_seq(self):
        sql = ('SELECT p.protein_acc, ps.sequence, pp.psm_id, peps.sequence '
               'FROM proteins AS p '
               'JOIN protein_seq AS ps USING(protein_acc) '
               'JOIN protein_psm AS pp USING(protein_acc) '
               'JOIN psms AS psms USING(psm_id) '
               'JOIN peptide_sequences AS peps USING(pep_id)'
               )
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_protein_psm_records(self):
        sql = ('SELECT protein_acc, psm_id FROM protein_psm')
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def get_protein_group_candidates(self):
        sql = ('SELECT pgm.master_id, pgm.psm_id, pp.protein_acc, '
               'peps.sequence, p.score, pev.evidence_lvl, pc.coverage '
               'FROM psm_protein_groups AS pgm '
               'JOIN protein_psm AS pp USING(psm_id) '
               'JOIN psms AS p USING(psm_id) '
               'JOIN peptide_sequences AS peps USING(pep_id) '
               'LEFT OUTER JOIN protein_evidence AS pev '
               'ON pev.protein_acc=pp.protein_acc '
               'LEFT OUTER JOIN protein_coverage AS pc '
               'ON pc.protein_acc=pp.protein_acc '
               'ORDER BY pgm.master_id'
               )
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def check_evidence_tables(self):
        """Returns True if there are records in evidence
        tables, otherwise returns False"""
        return self.check_table('protein_evidence') is not False

    def check_table(self, tablename):
        cursor = self.get_cursor()
        cursor.execute('SELECT * FROM {} LIMIT 10'.format(tablename))
        if len(cursor.fetchall()) > 0:
            return True
        return False

    def get_all_psms_proteingroups(self, evidence):
        fields = """
        SELECT pr.rownr, p.protein_acc, pgc.protein_acc, pgc.peptide_count,
        pgc.psm_count, pgc.protein_score, pc.coverage"""
        sql = """
        FROM psmrows AS pr
        JOIN psm_protein_groups AS ppg USING(psm_id)
        JOIN protein_group_master AS pgm USING(master_id)
        JOIN proteins AS p USING(pacc_id)
        JOIN protein_group_content AS pgc USING(master_id)
        JOIN protein_coverage AS pc ON pgc.protein_acc=pc.protein_acc
        """
        if evidence:
            sql = '{} {}'.format(sql, 
                    'JOIN protein_evidence AS pev ON pgc.protein_acc=pev.protein_acc')
            fields = '{}, pev.evidence_lvl '.format(fields)
        sql = '{} {} ORDER BY pr.rownr'.format(fields, sql)
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def select_all_psm_quants(self, shiftrows, isobaric=False, precursor=False):
        selects = ['pr.rownr']
        joins = ['JOIN psms USING(psm_id)', 'JOIN mzml USING(spectra_id)']
        sqlfields, fieldcount = {}, 1
        if isobaric:
            selects.extend(['ic.channel_name', 'iq.intensity', 'pif.pif'])
            joins.extend(['JOIN isobaric_quant AS iq USING(spectra_id)',
                          'JOIN isobaric_channels AS ic USING(channel_id)',
                          'LEFT OUTER JOIN precursor_ion_fraction AS pif USING(spectra_id)'])
            sqlfields['isochan'] = fieldcount
            sqlfields['isoquant'] = fieldcount + 1
            sqlfields['pif'] = fieldcount + 2
            fieldcount += 3
        if precursor:
            selects.extend(['pq.intensity', 'pfw.fwhm'])
            joins.extend(['LEFT OUTER JOIN ms1_align USING(spectra_id)',
                          'LEFT OUTER JOIN ms1_quant AS pq USING(feature_id)',
                          'LEFT OUTER JOIN ms1_fwhm AS pfw USING(feature_id)',
                          ])
            sqlfields['precursor'] = fieldcount
            sqlfields['fwhm'] = fieldcount + 1

        rowlim = 'WHERE pr.rownr>{} '.format(shiftrows - 1) if shiftrows else ''
        sql = ('SELECT {} FROM psmrows as pr {} {} '
               'ORDER BY pr.rownr'.format(', '.join(selects), ' '.join(joins), rowlim))
        cursor = self.get_cursor()
        return cursor.execute(sql), sqlfields

    def get_all_quantmaps(self):
        """Returns all unique quant channels from lookup as list"""
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT channel_name FROM isobaric_channels')
        return cursor.fetchall()

    def delete_sample_set_shift_rows(self, setnames):
        cursor = self.get_cursor()
        cursor.executemany('DELETE FROM biosets WHERE set_name=?', ((x,) for x in setnames))
        # Now all the rows will be gone where this set was, so we re-number:
        # Vacuuming updates the internal rowid column of the tables, which is
        # a count, when they are not INTEGER PRIMARY KEY
        self.conn.commit()
        self.conn.execute('DROP INDEX psmrowid_index')
        self.conn.execute('DROP INDEX psmrow_index')
        self.conn.execute('VACUUM')
        self.conn.execute('UPDATE psmrows SET rownr=rowid-1') # -1 since we start at 0
        self.conn.commit()
        # Put index back
        self.index_column('psmrowid_index', 'psmrows', 'psm_id')
        self.index_column('psmrow_index', 'psmrows', 'rownr')
        cursor = self.get_cursor()
        cursor.execute('DELETE FROM peptide_sequences WHERE pep_id IN ('
            'SELECT ps.pep_id FROM peptide_sequences AS ps '
            'LEFT OUTER JOIN psms ON ps.pep_id=psms.pep_id WHERE psms.pep_id IS NULL)')
