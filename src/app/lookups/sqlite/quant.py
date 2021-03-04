from app.lookups.sqlite.base import ResultLookupInterface


class QuantDB(ResultLookupInterface):

    def add_tables(self, tabletypes):
        if 'isobaric' in tabletypes:
            self.create_tables(['isobaric_quant', 'isobaric_channels',
                'precursor_ion_fraction'])
        if 'ms1' in tabletypes:
            self.create_tables(['ms1_quant', 'ms1_align', 'ms1_fwhm'])

    def get_fnfeats(self, fn_id):
        cursor = self.get_cursor()
        return cursor.execute(
            'SELECT mz, feature_id, charge, retention_time '
            'FROM ms1_quant '
            'WHERE mzmlfile_id=? ORDER BY mz', (fn_id,))

    def store_channelmap(self, channels):
        self.store_many(
            'INSERT OR IGNORE INTO isobaric_channels(channel_name) VALUES(?)', channels)

    def store_isobaric_quants(self, quants, pifs):
        self.store_many(
            'INSERT INTO isobaric_quant(spectra_id, channel_id, intensity) '
            'VALUES (?, ?, ?)', quants)
        self.store_many(
                'INSERT INTO precursor_ion_fraction(spectra_id, pif) VALUES (?, ?)', pifs)

    def index_isobaric_quants(self):
        self.index_column('spectraid_index', 'isobaric_quant', 'spectra_id')
        self.index_column('channel_id_index', 'isobaric_quant', 'channel_id')
        self.index_column('pif_spectraid_ix', 'precursor_ion_fraction', 'spectra_id')

    def get_specmap(self, fn_id, retention_time=False, scan_nr=False):
        """Returns all spectra ids for spectra filename, keyed by 
        retention time"""
        cursor = self.get_cursor()
        values = [fn_id]
        if retention_time:
            sql = 'SELECT retention_time,spectra_id FROM mzml WHERE mzmlfile_id=? '
        elif scan_nr:
            sql = 'SELECT scan_nr,spectra_id FROM mzml WHERE mzmlfile_id=? '
        cursor.execute(sql, tuple(values))
        return {k: {'id': sid} for k, sid in cursor.fetchall()}

    def get_channelmap(self):
        cursor = self.get_cursor()
        cursor.execute('SELECT channel_id, channel_name FROM isobaric_channels')
        return cursor

    def store_ms1_quants(self, quants):
        # Return feat IDs?
        quants = self.store_many_return_id(
            'INSERT INTO ms1_quant(mzmlfile_id, retention_time, mz, '
            'charge, intensity) VALUES (?, ?, ?, ?, ?)', quants)
        return quants

    def store_fwhm(self, quants):
        # FIXME get feature_id for all passed quants, possibly directly when storing them
        self.store_many('INSERT INTO ms1_fwhm(feature_id, fwhm) VALUES (?, ?)', quants)

    def store_ms1_alignments(self, aligns):
        self.store_many(
            'INSERT INTO ms1_align(spectra_id, feature_id) '
            'VALUES (?, ?)', aligns)

    def index_precursor_quants(self):
        self.index_column('ms1_feat_ix', 'ms1_quant', 'feature_id')
        self.index_column('charge_index', 'ms1_quant', 'charge')
        self.index_column('feat_mz_index', 'ms1_quant', 'mz')
        self.index_column('ms1_mzfn_index', 'ms1_quant', 'mzmlfile_id')
        self.index_column('fwhm_feat_index', 'ms1_fwhm', 'feature_id')

    def index_aligned_quants(self):
        self.index_column('ms1al_feat_ix', 'ms1_align', 'feature_id')
        self.index_column('ms1al_spec_ix', 'ms1_align', 'spectra_id')

    def get_spectra_mz_sorted(self):
        return self.get_cursor().execute(
            'SELECT spectra_id, mzmlfile_id, charge, mz, retention_time '
            'FROM mzml ORDER BY mzmlfile_id,mz')
