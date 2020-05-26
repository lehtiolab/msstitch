from app.lookups.sqlite.base import ResultLookupInterface


class ProtPepTable(ResultLookupInterface):
    table_map = {'protein': {'fntable': 'protein_tables',
                             'feattable': 'protein_group_master',
                             'isoqtable': 'protein_iso_quanted',
                             'isochtable': 'protquant_channels',
                             'prectable': 'protein_precur_quanted',
                             'fdrtable': 'protein_fdr',
                             },
                 'gene': {'fntable': 'gene_tables',
                          'feattable': 'genes',
                          'isoqtable': 'gene_iso_quanted',
                          'isochtable': 'genequant_channels',
                          'prectable': 'gene_precur_quanted',
                          'fdrtable': 'gene_fdr',
                          },
                 'assoc': {'fntable': 'gene_tables',
                           'feattable': 'associated_ids',
                           'isoqtable': 'assoc_iso_quanted',
                           'isochtable': 'genequant_channels',
                           'prectable': 'assoc_precur_quanted',
                           'fdrtable': 'assoc_fdr',
                           },
                 'peptide': {'fntable': 'peptide_tables',
                             'feattable': 'peptide_sequences',
                             'isoqtable': 'peptide_iso_quanted',
                             'isochtable': 'pepquant_channels',
                             'prectable': 'peptide_precur_quanted',
                             'fdrtable': 'peptide_fdr',
                             }
                 }

    def get_all_poolnames(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT set_name, set_id FROM biosets ORDER BY set_name')
        return cursor

    def store_table_files(self, tables):
        self.store_many(
            'INSERT INTO {}(set_id, filename) VALUES(?, ?)'.format(
                self.table_map[self.datatype]['fntable']),
            tables)

    def store_quant_channels(self, quantchannels, psmnrcols):
        table = self.table_map[self.datatype]['isochtable']
        if psmnrcols:
            sql = ('INSERT INTO {}({}, channel_name, amount_psms_name) VALUES'
                   '(?, ?, ?)'.format(table, self.colmap[table][1]))
        else:
            sql = ('INSERT INTO {}({}, channel_name) VALUES'
                   '(?, ?)'.format(table, self.colmap[table][1]))
        self.store_many(sql, quantchannels)

    def get_quantchannel_map(self):
        outdict = {}
        table = self.table_map[self.datatype]['isochtable']
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT channel_id, {}, channel_name, amount_psms_name'
            ' FROM {}'.format(self.colmap[table][1], table))
        for channel_id, fnid, channel_name, amount_psms_name in cursor:
            try:
                outdict[fnid][channel_name] = (channel_id, amount_psms_name)
            except KeyError:
                outdict[fnid] = {channel_name: (channel_id, amount_psms_name)}
        return outdict

    def store_isobaric_quants(self, quants, psmnrcols):
        table = self.table_map[self.datatype]['isoqtable']
        if psmnrcols:
            sql = ('INSERT INTO {}({}, channel_id, quantvalue, amount_psms) ' 
                   'VALUES ' '(?, ?, ?, ?)'.format(table, self.colmap[table][1]))
        else:
            sql = ('INSERT INTO {}({}, channel_id, quantvalue) ' 
                   'VALUES ' '(?, ?, ?)'.format(table, self.colmap[table][1]))
        self.store_many(sql, quants)

    def get_tablefn_map(self):
        table = self.table_map[self.datatype]['fntable']
        cursor = self.get_cursor()
        cursor.execute('SELECT * FROM {} '.format(table))
        return {fn: table_id for (table_id, setid, fn) in cursor}

    def get_feature_map(self):
        if self.datatype == 'protein':
            table = 'proteins'
        else:
            table = self.table_map[self.datatype]['feattable']
        columns = self.colmap[table][0:2]
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

    def update_selects(self, selectmap, fields, fieldcount):
        selectmap.update({field: i + fieldcount
                          for i, field in enumerate(fields)})
        fieldcount = max(selectmap.values()) + 1
        return selectmap, fieldcount

    def get_proteingroup_content(self):
        cursor = self.get_cursor()
        sql = ('SELECT pgm.protein_acc, pgc.protein_acc FROM '
               'protein_group_master AS pgm '
               'JOIN protein_group_content AS pgc USING(master_id)')
        return cursor.execute(sql)

    def get_isoquant_headernames(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT channel_name, amount_psms_name '
            'FROM {}'.format(self.table_map[self.datatype]['isochtable']))
        return cursor
