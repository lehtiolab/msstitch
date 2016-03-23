from app.lookups.sqlite.base import ResultLookupInterface


class BioSetDB(ResultLookupInterface):
    def add_tables(self):
        self.create_tables(['biosets', 'mzmlfiles'])

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
