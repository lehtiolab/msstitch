from app.lookups.sqlite.base import ResultLookupInterface


class ProtTableDB(ResultLookupInterface):
    def add_tables(self):
        self.create_tables(['protein_quanted', 'protquant_channels', 
                            'protein_precur_quanted'])

    def store_quant_channels(self, quantchannels):
        self.store_many(
            'INSERT INTO protquant_channels(protquant_file, channel_name, '
            'amount_psms_name) VALUES (?, ?, ?)',
            quantchannels)

    def get_protein_data(self, protein_acc):
        fields = ['psm.psm_id', 'psm.sequence', 'pgc.master',
                  'pgc.protein_acc', 'pcov.coverage', 'pd.description']
        joins = [('psm_protein_groups', 'ppg', 'master'),
                 ('protein_coverage', 'pcov', 'protein_acc'),
                 ('prot_desc', 'pd', 'protein_acc'),
                 ('psms', 'psm', 'psm_id'),
                 ]
        join_sql = '\n'.join(['JOIN {0} AS {1} USING({2})'.format(
            j[0], j[1], j[2]) for j in joins])
        sql = ('SELECT {0} FROM protein_group_content AS pgc {1}'
               'WHERE protein_acc="{2}"'.format(', '.join(fields),
                                                join_sql,
                                                protein_acc))
        cursor = self.get_cursor()
        return cursor.execute(sql).fetchall()

    def get_quanted_proteins(self, precursor=False, isobaric=False):
        if precursor and isobaric:
            sql = ('SELECT pq.protein_acc, pc.channel_name, pc.protquant_file, '
                   'pc.amount_psms_name, pq.quantvalue, pq.amount_psms, '
                   'preq.protquant_file, preq.quantvalue '
                   'FROM protein_quanted AS pq '
                   'JOIN protquant_channels AS pc USING(channel_id) '
                   'JOIN protein_precur_quanted AS preq USING(protein_acc) '
                   )
        elif precursor:
            sql = ('SELECT pq.protein_acc, pq.protquant_file, pq.quantvalue '
                   'FROM protein_precur_quanted AS pq')
        elif isobaric:
            sql = ('SELECT pq.protein_acc, pc.channel_name, pc.protquant_file, '
                   'pc.amount_psms_name, pq.quantvalue, pq.amount_psms '
                   'FROM protein_quanted AS pq '
                   'JOIN protquant_channels AS pc USING(channel_id) ')
        cursor = self.get_cursor()
        cursor.execute(sql)
        return cursor

    def get_quantchannel_headerfields(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT protquant_file, channel_name, amount_psms_name '
            'FROM protquant_channels')
        return cursor

    def get_precursorquant_headerfields(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT protquant_file '
            'FROM protein_precur_quanted')
        return cursor

    def get_quantchannel_map(self):
        outdict = {}
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT channel_id, protquant_file, channel_name, amount_psms_name'
            ' FROM protquant_channels')
        for channel_id, fname, channel_name, amount_psms_name in cursor:
            try:
                outdict[fname][channel_name] = (channel_id, amount_psms_name)
            except KeyError:
                outdict[fname] = {channel_name: (channel_id, amount_psms_name)}
        return outdict

    def store_isobaric_protquants(self, quants):
        self.store_many(
            'INSERT INTO protein_quanted(protein_acc, channel_id, quantvalue, '
            'amount_psms) VALUES (?, ?, ?, ?)', quants)

    def store_precursor_protquants(self, quants):
        self.store_many(
            'INSERT INTO protein_precur_quanted(protein_acc, protquant_file, quantvalue) '
            'VALUES (?, ?, ?)', quants)
