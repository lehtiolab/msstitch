from app.lookups.sqlite.base import DatabaseConnection


class SpectraDB(DatabaseConnection):
    def add_tables(self):
        self.create_tables(['mzml'])

    def store_mzmls(self, spectra):
        self.store_many(
            'INSERT INTO mzml(mzmlfilename, scan_nr, retention_time) '
            'VALUES (?, ?, ?)', spectra)

    def index_mzml(self):
        self.index_column('mzmlfn_index', 'mzml', 'mzmlfilename')
        self.index_column('scan_index', 'mzml', 'scan_nr')
        self.index_column('rt_index', 'mzml', 'retention_time')
