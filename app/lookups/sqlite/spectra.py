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
        self.index_column('mzmlfnid_index', 'mzml', 'mzmlfile_id')
        self.index_column('scan_index', 'mzml', 'scan_nr')
        self.index_column('rt_index', 'mzml', 'retention_time')
