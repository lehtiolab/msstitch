from app.lookups.sqlite.base import ResultLookupInterface
    

class QuantDB(ResultLookupInterface):

    def select_all_psm_quants(self):
        sql = (
    'SELECT pr.rownr, iq.quantmap, iq.intensity, pq.intensity '
    'FROM psmrows AS pr '
    'JOIN psms USING(psm_id) '
    'JOIN mzml USING(spectra_id) '
    'JOIN isobaric_quant AS iq USING(spectra_id) '
    'LEFT OUTER JOIN ms1_align USING(spectra_id) '
    'LEFT OUTER JOIN ms1_quant AS pq USING(feature_id)'
               )
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def lookup_retention_time(self, mzmlfn_id, scannr):
       # deprecate FIXME
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT retention_time FROM mzml '
            'WHERE mzmlfile_id=? AND scan_nr=?',
            (mzmlfn_id, scannr))
        return cursor.fetchall()

    def lookup_isobaric_quant(self, mzmlfn_id, scannr):
        # FIXME deprecate
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT iq.quantmap, iq.intensity '
            'FROM mzml AS mz '
            'JOIN isobaric_quant AS iq USING(spectra_id) '
            'WHERE mz.mzmlfile_id=? AND mz.scan_nr=?', (mzmlfn_id, scannr))
        return cursor.fetchall()

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
            'SELECT DISTINCT quantmap FROM isobaric_quant')
        return cursor.fetchall()


class IsobaricQuantDB(QuantDB):
    def add_tables(self):
        self.create_tables(['isobaric_quant'])

    def store_isobaric_quants(self, quants):
        self.store_many(
            'INSERT INTO isobaric_quant(spectra_id, quantmap, intensity) '
            'VALUES (?, ?, ?)', quants)

    def index_isobaric_quants(self):
        self.index_column('spectraid_index', 'isobaric_quant', 'spectra_id')


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
