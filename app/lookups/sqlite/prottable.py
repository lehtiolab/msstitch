from app.lookups.sqlite.base import ResultLookupInterface


class ProtPepTable(ResultLookupInterface):
    table_map = {'protein': {'fntable': 'protein_tables',
                             'feattable': 'proteins',
                             },
                 'peptide': {'fntable': 'peptide_tables',
                             'feattable': 'peptide_sequences'
                             }
                 }

    def store_table_files(self, tables):
        self.store_many(
            'INSERT INTO {}(set_id, filename) VALUES(?, ?)'.format(
                self.table_map[self.datatype]['fntable']),
            tables)

    def get_tablefn_map(self):
        table = self.table_map[self.datatype]['fntable']
        cursor = self.get_cursor()
        cursor.execute('SELECT * FROM {} '.format(table))
        return {fn: table_id for (table_id, setid, fn) in cursor}

    def get_feature_map(self):
        columns = {'protein': ['pacc_id', 'protein_acc'],
                   'peptide': ['pepid', 'sequence']
                   }
        table = self.table_map[self.datatype]['feattable']
        columns = columns[self.datatype]
        cursor = self.get_cursor()
        cursor.execute('SELECT {}, {} FROM {}'.format(columns[0], columns[1],
                                                      table))
        return {acc: table_id for (table_id, acc) in cursor}


class ProtTableDB(ProtPepTable):
    datatype = 'protein'

    def add_tables(self):
        self.create_tables(['protein_tables', 'protein_iso_quanted',
                            'protquant_channels', 'protein_precur_quanted',
                            'protein_probability', 'protein_fdr',
                            'protein_pep'])

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
                               probability=False, fdr=False, pep=False):
        selects = ['p.protein_acc']
        selectmap = {'p_acc': 0}
        selectfieldcount = max(selectmap.values()) + 1
        joins = []
        if isobaric:
            selects.extend(['pc.channel_name', 'bs.set_name',
                            'pc.amount_psms_name', 'piq.quantvalue',
                            'piq.amount_psms'])
            joins.extend([('protein_iso_quanted', 'piq', 'p', 'pacc_id', True),
                          ('protquant_channels', 'pc', 'piq', 'channel_id'),
                          ('protein_tables', 'pt', 'pc', 'prottable_id'),
                          ('biosets', 'bs', 'pt', 'set_id')
                          ])
            selectmap.update({field: i + selectfieldcount
                              for i, field in enumerate(
                                  ['channel', 'isoq_poolname',
                                   'isoq_psmsfield', 'isoq_val',
                                   'isoq_psms'])})
            selectfieldcount = max(selectmap.values()) + 1
        if precursor:
            selects.extend(['prqbs.set_name', 'preq.quantvalue'])
            joins.extend([('protein_precur_quanted', 'preq', 'p', 'pacc_id',
                           True),
                          ('protein_tables', 'prqpt', 'preq', 'prottable_id'),
                          ('biosets', 'prqbs', 'prqpt', 'set_id')
                          ])
            selectmap.update({field: i + selectfieldcount
                              for i, field in enumerate(['preq_poolname',
                                                         'preq_val'])})
            selectfieldcount = max(selectmap.values()) + 1
        if probability:
            selects.extend(['probbs.set_name', 'pprob.probability'])
            joins.extend([('protein_probability', 'pprob', 'p', 'pacc_id',
                           True),
                          ('protein_tables', 'probpt', 'pprob',
                           'prottable_id'),
                          ('biosets', 'probbs', 'probpt', 'set_id')
                          ])
            selectmap.update({field: i + selectfieldcount
                              for i, field in enumerate(['prob_poolname',
                                                         'prob_val'])})
            selectfieldcount = max(selectmap.values()) + 1
        if fdr:
            selects.extend(['fdrbs.set_name', 'pfdr.fdr'])
            joins.extend([('protein_fdr', 'pfdr', 'p', 'pacc_id', True),
                          ('protein_tables', 'fdrpt', 'pfdr', 'prottable_id'),
                          ('biosets', 'fdrbs', 'fdrpt', 'set_id')
                          ])
            selectmap.update({field: i + selectfieldcount
                              for i, field in enumerate(['fdr_poolname',
                                                         'fdr_val'])})
            selectfieldcount = max(selectmap.values()) + 1
        if pep:
            selects.extend(['pepbs.set_name', 'ppep.pep'])
            joins.extend([('protein_pep', 'ppep', 'p', 'pacc_id', True),
                          ('protein_tables', 'peppt', 'ppep', 'prottable_id'),
                          ('biosets', 'pepbs', 'peppt', 'set_id')
                          ])
            selectmap.update({field: i + selectfieldcount
                              for i, field in enumerate(['pep_poolname',
                                                         'pep_val'])})
            selectfieldcount = max(selectmap.values()) + 1

        sql = 'SELECT {} FROM proteins AS p'.format(
            ', '.join(selects))
        # NB Use full outer joins or left outer joins here on the stuff you
        # join to the proteins AS p table
        if joins:
            joinsql = ''
            for j in joins:
                joincmd = 'JOIN'
                if True in joins:
                    joincmd = 'LEFT OUTER {}'.format(joincmd)
                joinsql = '{5} {0} {1} AS {2} ON {3}.{4}={2}.{4}'.format(
                    joincmd, j[0], j[1], j[2], j[3], joinsql)
            sql = '{} {}'.format(sql, joinsql)
        sql = '{0} ORDER BY p.protein_acc'.format(sql)
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
            'INSERT INTO protein_iso_quanted(pacc_id, channel_id, '
            'quantvalue, amount_psms) '
            'VALUES (?, ?, ?, ?)', quants)

    def store_precursor_protquants(self, quants):
        self.store_many(
            'INSERT INTO protein_precur_quanted(pacc_id, prottable_id, '
            'quantvalue) '
            'VALUES (?, ?, ?)', quants)

    def store_protprob(self, probabilities):
        self.store_many(
            'INSERT INTO protein_probability(pacc_id, prottable_id, '
            'probability) '
            'VALUES (?, ?, ?)', probabilities)

    def store_protfdr(self, fdr):
        self.store_many(
            'INSERT INTO protein_fdr(pacc_id, prottable_id, '
            'fdr) '
            'VALUES (?, ?, ?)', fdr)

    def store_protpep(self, pep):
        self.store_many(
            'INSERT INTO protein_pep(pacc_id, prottable_id, '
            'pep) '
            'VALUES (?, ?, ?)', pep)
