from app.lookups.sqlite.base import ResultLookupInterface


class QuantDB(ResultLookupInterface):

    def select_all_psm_quants(self, isobaric=False, precursor=False):
        selects = ['pr.rownr']
        joins = ['JOIN psms USING(psm_id)', 'JOIN mzml USING(spectra_id)']
        sqlfields, fieldcount = {}, 1
        if isobaric:
            selects.extend(['ic.channel_name', 'iq.intensity'])
            joins.extend(['JOIN isobaric_quant AS iq USING(spectra_id)',
                          'JOIN isobaric_channels AS ic USING(channel_id)'])
            sqlfields['isochan'] = fieldcount
            sqlfields['isoquant'] = fieldcount + 1
            fieldcount += 2
        if precursor:
            selects.extend(['pq.intensity'])
            joins.extend(['LEFT OUTER JOIN ms1_align USING(spectra_id)',
                          'LEFT OUTER JOIN ms1_quant AS pq USING(feature_id)'])
            sqlfields['precursor'] = fieldcount
        if not precursor and not isobaric:
            raise RuntimeError('Cannot add quantification data, neither '
                               'isobaric, nor precursor have been specified.')
        sql = ('SELECT {} FROM psmrows as pr {} '
               'ORDER BY pr.rownr'.format(', '.join(selects), ' '.join(joins)))
        cursor = self.get_cursor()
        return cursor.execute(sql), sqlfields

    def get_precursor_quant_window(self, windowsize, minmz):
        cursor = self.get_cursor()
        return cursor.execute(
            'SELECT feature_id, mzmlfile_id, charge, mz, retention_time '
            'FROM ms1_quant '
            'WHERE mz > ? ORDER BY mz '
            'LIMIT ?', (minmz, windowsize))

    def get_all_quantmaps(self):
        """Returns all unique quant channels from lookup as list"""
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT channel_name FROM isobaric_channels')
        return cursor.fetchall()


class IsobaricQuantDB(QuantDB):
    def add_tables(self):
        self.create_tables(['isobaric_quant', 'isobaric_channels'])

    def store_channelmap(self, channels):
        self.store_many(
            'INSERT INTO isobaric_channels(channel_name) VALUES(?)', channels)

    def store_isobaric_quants(self, quants):
        self.store_many(
            'INSERT INTO isobaric_quant(spectra_id, channel_id, intensity) '
            'VALUES (?, ?, ?)', quants)

    def index_isobaric_quants(self):
        self.index_column('spectraid_index', 'isobaric_quant', 'spectra_id')
        self.index_column('channel_id_index', 'isobaric_quant', 'channel_id')

    def get_channelmap(self):
        cursor = self.get_cursor()
        cursor.execute('SELECT channel_id, channel_name FROM isobaric_channels')
        return cursor


class PrecursorQuantDB(QuantDB):
    def add_tables(self):
        self.create_tables(['ms1_quant', 'ms1_align'])

    def store_ms1_quants(self, quants):
        self.store_many(
            'INSERT INTO ms1_quant(mzmlfile_id, retention_time, mz, '
            'charge, intensity) VALUES (?, ?, ?, ?, ?)', quants)

    def store_ms1_alignments(self, aligns):
        self.store_many(
            'INSERT INTO ms1_align(spectra_id, feature_id) '
            'VALUES (?, ?)', aligns)

    def index_precursor_quants(self):
        self.index_column('charge_index', 'ms1_quant', 'charge')
        self.index_column('feat_mz_index', 'ms1_quant', 'mz')

    def get_spectra_mz_sorted(self):
        return self.get_cursor().execute(
            'SELECT spectra_id, mzmlfile_id, charge, mz, retention_time '
            'FROM mzml ORDER BY mz')
