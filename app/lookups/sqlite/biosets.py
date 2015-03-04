from app.lookups.sqlite.base import DatabaseConnection


class BioSetDB(DatabaseConnection):
    def add_tables(self):
        self.create_tables(['biosets', 'mzmlfiles'])

    def store_biosets(self, biosets):
        self.store_many('INSERT INTO biosets(set_name) VALUES (?)', biosets)

    def index_precursor_quants(self):
        self.index_column('setname_index', 'biosets', 'set_name')
