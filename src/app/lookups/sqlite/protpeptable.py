from app.lookups.sqlite.base import ResultLookupInterface


class ProtPepTable(ResultLookupInterface):
    table_map = {'protein': {'fntable': 'protein_tables',
                             'feattable': 'proteins',
                             'isoqtable': 'protein_iso_quanted',
                             'isochtable': 'protquant_channels',
                             'prectable': 'protein_precur_quanted',
                             'fdrtable': 'protein_fdr',
                             'fullqpsmtable': 'protein_iso_fullpsms',
                             },
                 'gene': {'fntable': 'gene_tables',
                          'feattable': 'genes',
                          'isoqtable': 'gene_iso_quanted',
                          'isochtable': 'genequant_channels',
                          'prectable': 'gene_precur_quanted',
                          'fdrtable': 'gene_fdr',
                          'fullqpsmtable': 'gene_iso_fullpsms',
                          },
                 'assoc': {'fntable': 'gene_tables',
                           'feattable': 'associated_ids',
                           'isoqtable': 'assoc_iso_quanted',
                           'isochtable': 'genequant_channels',
                           'prectable': 'assoc_precur_quanted',
                           'fdrtable': 'assoc_fdr',
                           'fullqpsmtable': 'assoc_iso_fullpsms',
                           },
                 'peptide': {'fntable': 'peptide_tables',
                             'feattable': 'peptide_sequences',
                             'isoqtable': 'peptide_iso_quanted',
                             'isochtable': 'pepquant_channels',
                             'prectable': 'peptide_precur_quanted',
                             'fdrtable': 'peptide_fdr',
                             'fullqpsmtable': 'peptide_iso_fullpsms',
                             'flrtable': 'ptm_flr',
                             }
                 }

    def __init__(self, fn=None):
        super().__init__(fn)
        self.singlecols_to_index = []

    def add_tables(self, tabletypes=[]):
        ttypes = ['fntable', 'isoqtable', 'isochtable', 'prectable', 'fdrtable',
                'fullqpsmtable', 'flrtable']
        self.create_tables([self.table_map[self.datatype][x] for x in ttypes
            if x in self.table_map[self.datatype]])

    def get_all_poolnames(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT set_name, set_id FROM biosets ORDER BY set_name')
        return cursor

    def store_table_files(self, tables):
        table = self.table_map[self.datatype]['fntable']
        self.store_many(
                'INSERT INTO {}(set_id, filename) VALUES(?, ?)'.format(table), tables)
        self.index_column('table_set_ix', table, 'set_id')

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

    def index_iso(self):
        table = self.table_map[self.datatype]['isoqtable']
        self.index_column('isoq_acc_ix', table, self.colmap[table][1])
        self.index_column('isoq_ch_ix', table, self.colmap[table][2])

    def index_singlecol_data(self):
        for ixname, table, col in self.singlecols_to_index:
            self.index_column(ixname, table, col)

    def get_tablefn_map(self):
        table = self.table_map[self.datatype]['fntable']
        cursor = self.get_cursor()
        cursor.execute('SELECT * FROM {} '.format(table))
        return {fn: table_id for (table_id, setid, fn) in cursor}

    def get_feature_map(self):
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
        table = self.table_map[self.datatype]['prectable']
        self.singlecols_to_index.append(('precursor_acc_ix', table, self.colmap[table][0]))
        self.singlecols_to_index.append(('precursor_table_ix', table, self.colmap[table][1]))

    def store_fdr(self, fdr):
        self.store_singlecol('fdrtable', fdr)
        table = self.table_map[self.datatype]['fdrtable']
        self.singlecols_to_index.append(('fdr_acc_ix', table, self.colmap[table][0]))
        self.singlecols_to_index.append(('fdr_table_ix', table, self.colmap[table][1]))

    def store_fullq_psms(self, fullqpsms):
        self.store_singlecol('fullqpsmtable', fullqpsms)
        table = self.table_map[self.datatype]['fullqpsmtable']
        self.singlecols_to_index.append(('fullq_acc_ix', table, self.colmap[table][0]))
        self.singlecols_to_index.append(('fullq_table_ix', table, self.colmap[table][1]))

    def store_ptm_flr(self, flrs):
        self.store_singlecol('flrtable', flrs)
        table = self.table_map[self.datatype]['flrtable']
        self.singlecols_to_index.append(('flr_acc_ix', table, self.colmap[table][0]))
        self.singlecols_to_index.append(('flr_table_ix', table, self.colmap[table][1]))

    def check_isoquant_psmnrs(self):
        cursor = self.get_cursor()
        tables = self.table_map[self.datatype]
        psmschecksql = """SELECT EXISTS(SELECT amount_psms FROM {} 
                WHERE amount_psms IS NOT NULL)""".format(tables['isoqtable'])
        return cursor.execute(psmschecksql).fetchone()[0]

    def get_isoquant_headernames(self, stored_psmnrs):
        cursor = self.get_cursor()
        tables = self.table_map[self.datatype]
        if stored_psmnrs:
            cursor.execute(
                'SELECT DISTINCT channel_name, amount_psms_name '
                'FROM {}'.format(tables['isochtable']))
        else:
            cursor.execute(
                'SELECT DISTINCT channel_name FROM {}'.format(tables['isochtable']))
        return cursor
