from app.lookups.sqlite.biosets import BioSetDB


class SpectraDB(BioSetDB):
    def add_tables(self):
        super().add_tables()
        self.create_tables(['mzml'])

    def store_mzmls(self, spectra):
        self.store_many(
            'INSERT INTO mzml(mzmlfile_id, scan_nr, retention_time) '
            'VALUES (?, ?, ?)', spectra)

    def index_mzml(self):
        self.index_column('spectra_id_index', 'mzml', 'spectra_id')
        self.index_column('mzmlfnid_mzml_index', 'mzml', 'mzmlfile_id')
        self.index_column('scan_index', 'mzml', 'scan_nr')
        self.index_column('rt_index', 'mzml', 'retention_time')

    def get_spectradata(self, scannr):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT sp.mzmlfile_id, bs.set_name, sp.retention_time '
            'FROM mzml AS sp '
            'JOIN mzmlfiles AS mf USING(mzmlfile_id) '
            'JOIN biosets AS bs USING(set_id) '
            'WHERE sp.scan_nr=?'
            , (scannr,))
        return cursor
