from app.lookups.sqlite.base import DatabaseConnection


class QuantDB(DatabaseConnection):
    def create_quantdb(self, workdir):
        self.create_db(workdir,
                       {'quant': ['spectra_filename TEXT', 'scan_nr TEXT',
                                  'quantmap TEXT', 'intensity REAL']})

    def store_quants(self, quants):
        cursor = self.get_cursor()
        cursor.executemany(
            'INSERT INTO quant(spectra_filename, scan_nr, quantmap, intensity)'
            ' VALUES (?, ?, ?, ?)',
            quants)
        self.conn.commit()

    def index_quants(self):
        self.index_column('spec_index', 'quant', 'spectra_filename')
        self.index_column('scan_index', 'quant', 'scan_nr')

    def lookup_quant(self, spectrafile, scannr):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT quantmap, intensity FROM quant WHERE '
            'spectra_filename=? AND scan_nr=?', (spectrafile, scannr))
        return cursor.fetchall()

    def get_all_quantmaps(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT quantmap FROM quant')
        return cursor.fetchall()
