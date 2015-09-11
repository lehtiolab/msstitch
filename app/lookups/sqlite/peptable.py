from app.lookups.sqlite.protpeptable import ProtPepTable


class PepTableDB(ProtPepTable):
    datatype = 'peptide'
    colmap = {'peptide_precur_quanted': ['pep_id', 'peptable_id', 'quant'],
              'peptide_fdr': ['pep_id', 'peptable_id', 'fdr'],
              'peptide_pep': ['pep_id', 'peptable_id', 'pep'],
              }

    def add_tables(self):
        self.create_tables(['peptide_tables', 'pepquant_channels',
                            'peptide_iso_quanted', 'peptide_precur_quanted',
                            'peptide_fdr', 'peptide_pep'])

    def get_isoquant_channels(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT channel_name '
            'FROM pepquant_channels')
        return (x[0] for x in cursor)

    def store_quant_channels(self, quantchannels):
        self.store_many(
            'INSERT INTO pepquant_channels(peptable_id, channel_name) '
            'VALUES (?, ?)',
            quantchannels)

    def store_isobaric_quants(self, quants):
        self.store_many(
            'INSERT INTO peptide_iso_quanted(pep_id, channel_id, quantvalue) '
            'VALUES (?, ?, ?)', quants)

    def get_quantchannel_map(self):
        outdict = {}
        amount_psms_name = None
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT channel_id, peptable_id, channel_name '
            'FROM pepquant_channels')
        for channel_id, fnid, channel_name in cursor:
            try:
                outdict[fnid][channel_name] = (channel_id, amount_psms_name)
            except KeyError:
                outdict[fnid] = {channel_name: (channel_id, amount_psms_name)}
        return outdict

    def prepare_mergetable_sql(self, precursor=False, isobaric=False,
                               fdr=False, pep=False):
        selects = ['p.sequence']
        selectmap, count = self.update_selects({}, ['p_seq'], 0)
        joins = []
        if isobaric:
            selects.extend(['pc.channel_name', 'bs.set_name',
                            'piq.quantvalue'])
            joins.extend([('peptide_iso_quanted', 'piq', 'p', 'pep_id', True),
                          ('pepquant_channels', 'pc', 'piq', 'channel_id'),
                          ('peptide_tables', 'pt', 'pc', 'peptable_id'),
                          ('biosets', 'bs', 'pt', 'set_id')
                          ])
            fld = ['channel', 'isoq_poolname', 'isoq_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if precursor:
            selects.extend(['prqbs.set_name', 'preq.quant'])
            joins.extend([('peptide_precur_quanted', 'preq', 'p', 'pep_id',
                           True),
                          ('peptide_tables', 'prqpt', 'preq', 'peptable_id'),
                          ('biosets', 'prqbs', 'prqpt', 'set_id')
                          ])
            fld = ['preq_poolname', 'preq_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if fdr:
            selects.extend(['fdrbs.set_name', 'pfdr.fdr'])
            joins.extend([('peptide_fdr', 'pfdr', 'p', 'pep_id', True),
                          ('peptide_tables', 'fdrpt', 'pfdr', 'peptable_id'),
                          ('biosets', 'fdrbs', 'fdrpt', 'set_id')
                          ])
            fld = ['fdr_poolname', 'fdr_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if pep:
            selects.extend(['pepbs.set_name', 'ppep.pep'])
            joins.extend([('peptide_pep', 'ppep', 'p', 'pep_id', True),
                          ('peptide_tables', 'peppt', 'ppep', 'peptable_id'),
                          ('biosets', 'pepbs', 'peppt', 'set_id')
                          ])
            fld = ['pep_poolname', 'pep_val']
            selectmap, count = self.update_selects(selectmap, fld, count)

        sql = 'SELECT {} FROM peptide_sequences AS p'.format(
            ', '.join(selects))
        sql = self.get_sql_joins_mergetable(sql, joins)
        sql = '{0} ORDER BY p.sequence'.format(sql)
        return sql, selectmap
