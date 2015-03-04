from app.lookups.sqlite.base import DatabaseConnection


class BioSetDB(DatabaseConnection):
    def add_tables(self):
        self.create_tables(['biosets', 'mzmlfiles'])

    def store_biosets(self, biosets):
        self.store_many('INSERT INTO biosets(set_name) VALUES (?)', biosets)

    def store_mzmlfiles(self, mzmlfiles):
        self.store_many('INSERT INTO mzmlfiles(mzmlfilname, set_id) '
                        'VALUES (?)', mzmlfiles)

    def index_biosets(self):
        self.index_column('setname_index', 'biosets', 'set_name')
        self.index_column('mzmlfn_index', 'mzmlfiles', 'mzmlfilename')
