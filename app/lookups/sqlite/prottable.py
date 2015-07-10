from app.lookups.sqlite.base import ResultLookupInterface


class ProtTableDB(ResultLookupInterface):
    def add_tables(self):
        self.create_tables(['protein_tables', 'protein_quanted',
                            'protquant_channels', 'protein_precur_quanted',
                            'protein_probability'])

    def store_protein_tables(self, tables):
        self.store_many(
            'INSERT INTO protein_tables(set_id, prottable_file) VALUES(?, ?)',
            tables)

    def store_quant_channels(self, quantchannels):
        self.store_many(
            'INSERT INTO protquant_channels(prottable_id, channel_name, '
            'amount_psms_name) VALUES (?, ?, ?)',
            quantchannels)

    def get_all_protein_psms_with_sets(self):
        return self.get_proteins_psms(extended=True)

    def get_all_proteins_psms_for_unipeps(self):
        return self.get_proteins_psms(extended=False)

    def get_proteins_psms(self, extended=False):
        fields = ['pgm.protein_acc', 'sets.set_name',
                  'psm.sequence']
        joins = [('psm_protein_groups', 'ppg', 'master_id'),
                 ('psms', 'psm', 'psm_id'),
                 ('mzml', 'sp', 'spectra_id'),
                 ('mzmlfiles', 'mzfn', 'mzmlfile_id'),
                 ('biosets', 'sets', 'set_id'),
                 ]
        if extended:
            fields.extend(['psm.psm_id', 'pd.description', 'pcov.coverage'])
            joins.extend([('prot_desc', 'pd', 'protein_acc'),
                          ('protein_coverage', 'pcov', 'protein_acc'),
                          ])
        sql = ('SELECT {} FROM protein_group_master '
               'AS pgm'.format(', '.join(fields)))
        join_sql = ' '.join(['JOIN {} AS {} USING({})'.format(
            j[0], j[1], j[2]) for j in joins])
        sql = '{} {} ORDER BY pgm.protein_acc, sets.set_name'.format(sql,
                                                                     join_sql)
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def prepare_mergetable_sql(self, precursor=False, isobaric=False,
                               probability=False):
        selects = ['pq.protein_acc']
        selectmap = {'p_acc': 0}
        selectfieldcount = max(selectmap.values()) + 1
        joins = []
        if isobaric:
            selects.extend(['pc.channel_name', 'bs.set_name',
                            'pc.amount_psms_name', 'pq.quantvalue',
                            'pq.amount_psms'])
            joins.extend([('protquant_channels', 'pc', 'channel_id'),
                          ('protein_tables', 'pt', 'prottable_id'),
                          ('biosets', 'bs', 'set_id')
                          ])
            selectmap.update({field: i + selectfieldcount
                              for i, field in enumerate(
                                  ['channel', 'isoq_poolname',
                                   'isoq_psmsfield', 'isoq_val',
                                   'isoq_psms'])})
            selectfieldcount = max(selectmap.values()) + 1
        if precursor:
            selects.extend(['prqbs.set_name', 'preq.quantvalue'])
            joins.extend([('protein_precur_quanted', 'preq', 'protein_acc'),
                          ('protein_tables', 'prqpt', 'prottable_id'),
                          ('biosets', 'prqbs', 'set_id')
                          ])
            selectmap.update({field: i + selectfieldcount
                              for i, field in enumerate(['preq_poolname',
                                                         'preq_val'])})
            selectfieldcount = max(selectmap.values()) + 1
        if probability:
            selects.extend(['probbs.set_name', 'pprob.probability'])
            joins.extend([('protein_probability', 'pprob', 'protein_acc'),
                          ('protein_tables', 'probpt', 'prottable_id'),
                          ('biosets', 'probbs', 'set_id')
                          ])
            selectmap.update({field: i + selectfieldcount
                              for i, field in enumerate(['prob_poolname',
                                                         'prob_val'])})
            selectfieldcount = max(selectmap.values()) + 1

        sql = 'SELECT {} FROM protein_quanted AS pq'.format(', '.join(selects))
        if joins:
            sql = '{} {}'.format(sql, ' '.join(
                ['JOIN {} AS {} USING({})'.format(j[0], j[1], j[2])
                 for j in joins]))
        sql = '{0} ORDER BY pq.protein_acc'.format(sql)
        return sql, selectmap

    def get_merged_proteins(self, sql):
        cursor = self.get_cursor()
        cursor.execute(sql)
        return cursor

    def get_isoquant_amountpsms_channels(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT channel_name, amount_psms_name '
            'FROM protquant_channels')
        return cursor

    def get_precursorquant_headerfields(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT prottable_id '
            'FROM protein_precur_quanted')
        return cursor

    def get_all_poolnames(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT set_name, set_id '
            'FROM biosets')
        return cursor

    def get_protein_table_map(self):
        cursor = self.get_cursor()
        cursor.execute('SELECT prottable_id, prottable_file '
                       'FROM protein_tables')
        return {fn: table_id for (table_id, fn) in cursor}

    def get_quantchannel_map(self):
        outdict = {}
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT channel_id, prottable_id, channel_name, amount_psms_name'
            ' FROM protquant_channels')
        for channel_id, fnid, channel_name, amount_psms_name in cursor:
            try:
                outdict[fnid][channel_name] = (channel_id, amount_psms_name)
            except KeyError:
                outdict[fnid] = {channel_name: (channel_id, amount_psms_name)}
        return outdict

    def store_isobaric_protquants(self, quants):
        self.store_many(
            'INSERT INTO protein_quanted(protein_acc, channel_id, quantvalue, '
            'amount_psms) '
            'VALUES (?, ?, ?, ?)', quants)

    def store_precursor_protquants(self, quants):
        self.store_many(
            'INSERT INTO protein_precur_quanted(protein_acc, prottable_id, '
            'quantvalue) '
            'VALUES (?, ?, ?)', quants)

    def store_protprob(self, probabilities):
        self.store_many(
            'INSERT INTO protein_probability(protein_acc, prottable_id, '
            'probability) '
            'VALUES (?, ?, ?)', probabilities)
