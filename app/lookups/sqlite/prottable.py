from app.lookups.sqlite.base import ResultLookupInterface


class ProtPepTable(ResultLookupInterface):
    table_map = {'protein': {'fntable': 'protein_tables',
                             'feattable': 'proteins',
                             'prectable': 'protein_precur_quanted',
                             'fdrtable': 'protein_fdr',
                             'peptable': 'protein_pep',
                             },
                 'peptide': {'fntable': 'peptide_tables',
                             'feattable': 'peptide_sequences',
                             'prectable': 'peptide_precur_quanted',
                             'fdrtable': 'peptide_fdr',
                             'peptable': 'peptide_pep',
                             }
                 }

    def get_all_poolnames(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT set_name, set_id FROM biosets')
        return cursor

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
                   'peptide': ['pep_id', 'sequence']
                   }
        table = self.table_map[self.datatype]['feattable']
        columns = columns[self.datatype]
        cursor = self.get_cursor()
        cursor.execute('SELECT {}, {} FROM {}'.format(columns[0], columns[1],
                                                      table))
        return {acc: table_id for (table_id, acc) in cursor}

    def store_singlecol(self, tablekey, vals):
        table = self.table_map[self.datatype][tablekey]
        cols = self.colmap[table]
        self.store_many('INSERT INTO {}({}, {}, {}) '
                        'VALUES (?, ?, ?)'.format(table, cols[0], cols[1],
                                                  cols[2]), vals)

    def store_precursor_quants(self, quants):
        self.store_singlecol('prectable', quants)

    def store_fdr(self, fdr):
        self.store_singlecol('fdrtable', fdr)

    def store_pep(self, pep):
        self.store_singlecol('peptable', pep)


class PepTableDB(ProtPepTable):
    datatype = 'peptide'
    colmap = {'peptide_precur_quanted': ['pep_id', 'peptable_id', 'quant'],
              'peptide_fdr': ['pep_id', 'peptable_id', 'fdr'],
              'peptide_pep': ['pep_id', 'peptable_id', 'pep'],
              }

    def add_tables(self):
        self.create_tables(['peptide_tables', 'pepquant_channels',
                            'peptide_iso_quanted', 'peptide_precur_quanted',
                            'peptide_fdr', 'peptide_pep'])

    def store_quant_channels(self, quantchannels):
        self.store_many(
            'INSERT INTO pepquant_channels(peptable_id, channel_name) '
            'VALUES (?, ?)',
            quantchannels)

    def store_isobaric_quants(self, quants):
        self.store_many(
            'INSERT INTO peptide_iso_quanted(pacc_id, channel_id, quantvalue) '
            'VALUES (?, ?, ?)', quants)

    def get_quantchannel_map(self):
        outdict = {}
        amount_psms_name = None
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT channel_id, prottable_id, channel_name '
            'FROM protquant_channels')
        for channel_id, fnid, channel_name in cursor:
            try:
                outdict[fnid][channel_name] = (channel_id, amount_psms_name)
            except KeyError:
                outdict[fnid] = {channel_name: (channel_id, amount_psms_name)}
        return outdict


class ProtTableDB(ProtPepTable):
    datatype = 'protein'
    colmap = {'protein_precur_quanted': ['pacc_id', 'prottable_id', 'quant'],
              'protein_fdr': ['pacc_id', 'prottable_id', 'fdr'],
              'protein_pep': ['pacc_id', 'prottable_id', 'pep'],
              'protein_probability': ['pacc_id', 'prottable_id',
                                      'probability'],
              }

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
                  'pep.sequence']
        joins = [('psm_protein_groups', 'ppg', 'master_id'),
                 ('psms', 'psm', 'psm_id'),
                 ('peptide_sequences', 'pep', 'pep_id'),
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

    def update_selects(self, selectmap, fields, fieldcount):
        selectmap.update({field: i + fieldcount
                          for i, field in enumerate(fields)})
        fieldcount = max(selectmap.values()) + 1
        return selectmap, fieldcount

    def prepare_mergetable_sql(self, precursor=False, isobaric=False,
                               probability=False, fdr=False, pep=False):
        selects = ['p.protein_acc']
        selectmap, count = self.update_selects({}, ['p_acc'], 0)
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
            fld = ['channel', 'isoq_poolname', 'isoq_psmsfield', 'isoq_val',
                   'isoq_psms']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if precursor:
            selects.extend(['prqbs.set_name', 'preq.quant'])
            joins.extend([('protein_precur_quanted', 'preq', 'p', 'pacc_id',
                           True),
                          ('protein_tables', 'prqpt', 'preq', 'prottable_id'),
                          ('biosets', 'prqbs', 'prqpt', 'set_id')
                          ])
            fld = ['preq_poolname', 'preq_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if probability:
            selects.extend(['probbs.set_name', 'pprob.probability'])
            joins.extend([('protein_probability', 'pprob', 'p', 'pacc_id',
                           True),
                          ('protein_tables', 'probpt', 'pprob',
                           'prottable_id'),
                          ('biosets', 'probbs', 'probpt', 'set_id')
                          ])
            fld = ['prob_poolname', 'prob_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if fdr:
            selects.extend(['fdrbs.set_name', 'pfdr.fdr'])
            joins.extend([('protein_fdr', 'pfdr', 'p', 'pacc_id', True),
                          ('protein_tables', 'fdrpt', 'pfdr', 'prottable_id'),
                          ('biosets', 'fdrbs', 'fdrpt', 'set_id')
                          ])
            fld = ['fdr_poolname', 'fdr_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if pep:
            selects.extend(['pepbs.set_name', 'ppep.pep'])
            joins.extend([('protein_pep', 'ppep', 'p', 'pacc_id', True),
                          ('protein_tables', 'peppt', 'ppep', 'prottable_id'),
                          ('biosets', 'pepbs', 'peppt', 'set_id')
                          ])
            fld = ['pep_poolname', 'pep_val']
            selectmap, count = self.update_selects(selectmap, fld, count)

        sql = 'SELECT {} FROM proteins AS p'.format(
            ', '.join(selects))
        sql = self.get_sql_joins_mergetable(sql, joins)
        return sql, selectmap

    def get_sql_joins_mergetable(self, sql, joins):
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
        return sql

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

    def store_isobaric_quants(self, quants):
        self.store_many(
            'INSERT INTO protein_iso_quanted(pacc_id, channel_id, '
            'quantvalue, amount_psms) '
            'VALUES (?, ?, ?, ?)', quants)

    def store_protprob(self, probabilities):
        self.store_many(
            'INSERT INTO protein_probability(pacc_id, prottable_id, '
            'probability) '
            'VALUES (?, ?, ?)', probabilities)
