from app.lookups.sqlite.base import ResultLookupInterface


class ProtTableDB(ResultLookupInterface):
    def add_tables(self):
        self.create_tables(['protein_quanted', 'protquant_channels',
                            'protein_quanted_psms'])

    def store_quant_channels(self, quantchannels):
        self.store_many(
            'INSERT INTO protquant_channels(channel_name) VALUES (?)',
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
               'WHERE master IN '
               '(SELECT master FROM protein_group_content '
               'WHERE protein_acc="{2}")'.format(', '.join(fields),
                                                 join_sql,
                                                 protein_acc))
        cursor = self.get_cursor()
        return cursor.execute(sql).fetchall()

    def get_quantchannel_ids(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT channel_id, channel_name FROM protquant_channels')
        return {channel: chan_id for chan_id, channel in cursor}

    def store_protquants(self, quants):
        self.store_many(
            'INSERT INTO protein_quanted(protein_acc, channel_id, quantvalue) '
            'VALUES (?, ?, ?)', quants)

