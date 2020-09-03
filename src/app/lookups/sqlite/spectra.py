from app.lookups.sqlite.base import ResultLookupInterface


class SpectraDB(ResultLookupInterface):
    def add_tables(self, tabletypes):
        self.create_tables(['biosets', 'mzmlfiles'])
        self.create_tables(['mzml', 'ioninjtime', 'ionmob'])

    def store_biosets(self, biosets):
        self.store_many('INSERT INTO biosets(set_name) VALUES (?)', biosets)

    def store_mzmlfiles(self, mzmlfiles):
        self.store_many('INSERT INTO mzmlfiles(mzmlfilename, set_id) '
                        'VALUES (?, ?)', mzmlfiles)

    def index_biosets(self):
        self.index_column('setname_index', 'biosets', 'set_name')
        self.index_column('setid_index', 'biosets', 'set_id')
        self.index_column('setid_mzmlfn_index', 'mzmlfiles', 'set_id')
        self.index_column('mzmlfn_index', 'mzmlfiles', 'mzmlfilename')
        self.index_column('mzmlfnid_index', 'mzmlfiles', 'mzmlfile_id')

    def get_setnames(self):
        cursor = self.get_cursor()
        cursor.execute('SELECT set_name, set_id FROM biosets')
        return {setname: set_id for setname, set_id in cursor.fetchall()}
    def store_mzmls(self, spectra, ioninj, ionmob):
        self.store_many(
            'INSERT INTO mzml(spectra_id, mzmlfile_id, scan_sid, charge, mz, '
            'retention_time) '
            'VALUES (?, ?, ?, ?, ?, ?)', spectra)
        self.store_many(
            'INSERT INTO ioninjtime(spectra_id, ion_injection_time) '
            'VALUES(?, ?)', ioninj)
        self.store_many(
            'INSERT INTO ionmob(spectra_id, ion_mobility) '
            'VALUES(?, ?)', ionmob)

    def index_mzml(self):
        self.index_column('spectra_id_index', 'mzml', 'spectra_id')
        self.index_column('mzmlfnid_mzml_index', 'mzml', 'mzmlfile_id')
        self.index_column('scan_index', 'mzml', 'scan_sid')
        self.index_column('specrt_index', 'mzml', 'retention_time')
        self.index_column('specmz_index', 'mzml', 'mz')
        self.index_column('ionm_sp_ix', 'ionmob', 'spectra_id')
        self.index_column('ionin_sp_ix', 'ioninjtime', 'spectra_id')
