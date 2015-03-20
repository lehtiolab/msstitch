from app.lookups.sqlite.base import ResultLookupInterface


class ProtQuantDB(ResultLookupInterface):
    def add_tables(self):
        self.create_tables(['protein_quanted', 'protquant_channels',
                            'protein_quanted_psms'])

    def store_quant_channels(self, quantchannels):
        self.store_many(
            'INSERT INTO protquant_channels(channel_name) VALUES (?)',
            quantchannels)

    def get_quantchannel_ids(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT channel_id, channel_name FROM protquant_channels')
        return {channel: chan_id for chan_id, channel in cursor}

    def store_protquants(self, quants):
        self.store_many(
            'INSERT INTO protein_quanted(protein_acc, channel_id, quantvalue) '
            'VALUES (?, ?, ?)', quants)

