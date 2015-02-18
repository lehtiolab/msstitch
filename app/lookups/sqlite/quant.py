from app.lookups.sqlite.base import DatabaseConnection


class QuantDB(DatabaseConnection):
    def create_quantdb(self, workdir):
        self.create_db(workdir,
                       {
                        'mzml': ['mzmlfilename TEXT',
                                 'scan_nr TEXT',
                                 'retention_time REAL'],
                        # FIXME in future, make possible to build full db in
                        # small steps. This will speed up pipeline since
                        # identification takes longer time than quanting. Also
                        # would be good to make one large sqlite file with id,
                        # protein group, quant in it. Problem is: ID is
                        # slowest, but thats where the foreign keys will refer
                        # to, if we dont do a RT/MS2/file db from the mzml
                        # FIXME all iso/precursor quants, or only ones with IDs?
                        # if only ones with ids, we can foreign key to ID
                        # table.
                        # FIXME isoquant should get retention time instead of
                        # scan nr so we can fk join on rt/specfile to the
                        # mzml table
                        # FIXME lookup of ms1 quant like this:
                        # scannr -> rt from mzmltable; rt interval & mz
                        # interval -> ms1quants -> select best match on mz
                        'isobaric_quant': ['mzmlfilename TEXT',
                                           'retention_time REAL',
                                           'quantmap TEXT',
                                           'intensity REAL',
                                           'FOREIGN KEY(mzmlfilename)'
                                           'REFERENCES mzml '
                                           'FOREIGN KEY(retention_time)'
                                           'REFERENCES mzml '
                                           ],
                        'ms1_quant': ['mzmlfilename TEXT',
                                      'retention_time REAL', 'mz REAL',
                                      'charge INTEGER', 'intensity REAL',
                                      'FOREIGN KEY(mzmlfilename)'
                                      'REFERENCES mzml '
                                      'FOREIGN KEY(retention_time)'
                                      'REFERENCES mzml ']
                                      }, foreign_keys=True)

    def store_isobaric_quants(self, quants):
        self.store_many(
            'INSERT INTO isobaric_quant(mzmlfilename, retention_time, '
            'quantmap, intensity) VALUES (?, ?, ?, ?)', quants)

    def store_mzmls(self, spectra):
        self.store_many(
            'INSERT INTO mzml(mzmlfilename, scan_nr, retention_time) '
            'VALUES (?, ?, ?)', spectra)

    def store_ms1_quants(self, quants):
        self.store_many(
            'INSERT INTO ms1_quant(mzmlfilename, retention_time, mz, '
            'charge, intensity) VALUES (?, ?, ?, ?, ?)', quants)

    def store_many(self, sql, values):
        cursor = self.get_cursor()
        cursor.executemany(sql, values)
        self.conn.commit()

    def index_mzml(self):
        self.index_column('mzmlfn_index', 'mzml', 'mzmlfilename')
        self.index_column('scan_index', 'mzml', 'scan_nr')
        self.index_column('rt_index', 'mzml', 'retention_time')

    def index_isobaric_quants(self):
        pass

    def index_precursor_quants(self):
        self.index_column('charge_index', 'ms1_quant', 'charge')
        self.index_column('mz_index', 'ms1_quant', 'mz')

    def lookup_retention_time(self, spectrafile, scannr):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT retention_time FROM mzml '
            'WHERE mzmlfilename=? AND scan_nr=?',
            (spectrafile, scannr))
        return cursor.fetchall()

    def lookup_isobaric_quant(self, spectrafile, scannr):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT iq.quantmap, iq.intensity '
            'FROM mzml AS mz '
            'JOIN isobaric_quant AS iq USING(retention_time) '
            'WHERE mz.mzmlfilename=? AND mz.scan_nr=?', (spectrafile, scannr))
        return cursor.fetchall()

    def lookup_precursor_quant(self, spectrafile, charge, minrt, maxrt,
                               minmz, maxmz):
        # FIXME check if this is slow since it has two BETWEEN in it
        # we could replace the mz BETWEEN by filtering in python
        cursor = self.get_cursor()
        return cursor.execute(
            'SELECT retention_time, mz, intensity '
            'FROM precur_quants '
            'WHERE mzmlfilename=? AND charge=? AND '
            'rt BETWEEN ? AND ? AND mz BETWEEN ? AND ?',
            (spectrafile, charge, minrt, maxrt, minmz, maxmz))

    def get_all_quantmaps(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT quantmap FROM isobaric_quant')
        return cursor.fetchall()
