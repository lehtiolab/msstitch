from app.lookups.sqlite.protpeptable import ProtPepTable


class ProtTableDB(ProtPepTable):
    datatype = 'protein'
    colmap = {'protein_precur_quanted': ['pacc_id', 'prottable_id', 'quant'],
              'protein_fdr': ['pacc_id', 'prottable_id', 'fdr'],
              'protein_pep': ['pacc_id', 'prottable_id', 'pep'],
              'protein_probability': ['pacc_id', 'prottable_id',
                                      'probability'],
              }

    def add_tables(self):
        self.create_tables(['protein_tables', 'protein_iso_quanted',
                            'protquant_channels', 'protein_precur_quanted',
                            'protein_probability', 'protein_fdr',
                            'protein_pep'])

    def store_quant_channels(self, quantchannels):
        self.store_many(
            'INSERT INTO protquant_channels(prottable_id, channel_name, '
            'amount_psms_name) VALUES (?, ?, ?)',
            quantchannels)

    def get_all_protein_psms_with_sets(self):
        return self.get_proteins_psms(extended=True)

    def get_all_proteins_psms_for_unipeps(self):
        return self.get_proteins_psms(extended=False)

    def prepare_mergetable_sql(self, precursor=False, isobaric=False,
                               probability=False, fdr=False, pep=False):
        selects = ['p.protein_acc', 'bs.set_name']
        selectmap, count = self.update_selects({}, ['p_acc', 'set_name'], 0)
        joins = []
        if isobaric:
            selects.extend(['pc.channel_name',
                            'pc.amount_psms_name', 'piq.quantvalue',
                            'piq.amount_psms'])
            joins.extend([('protquant_channels', 'pc', ['pt']),
                          ('protein_iso_quanted', 'piq', ['p', 'pc'], True),
                          ])
            fld = ['channel', 'isoq_psmsfield', 'isoq_val',
                   'isoq_psms']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if precursor:
            selects.extend(['preq.quant'])
            joins.append(('protein_precur_quanted', 'preq', ['p', 'pt'], True))
            fld = ['preq_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if probability:
            selects.extend(['pprob.probability'])
            joins.append(('protein_probability', 'pprob', ['p', 'pt'], True))
            fld = ['prob_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if fdr:
            selects.extend(['pfdr.fdr'])
            joins.append(('protein_fdr', 'pfdr', ['p', 'pt'], True))
            fld = ['fdr_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if pep:
            selects.extend(['ppep.pep'])
            joins.append(('protein_pep', 'ppep', ['p', 'pt'], True))
            fld = ['pep_val']
            selectmap, count = self.update_selects(selectmap, fld, count)

        sql = ('SELECT {} FROM proteins AS p JOIN biosets AS bs '
               'JOIN protein_tables AS pt'.format(', '.join(selects)))
        sql = self.get_sql_joins_mergetable(sql, joins)
        sql = '{0} ORDER BY p.protein_acc'.format(sql)
        return sql, selectmap

    def get_isoquant_amountpsms_channels(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT channel_name, amount_psms_name '
            'FROM protquant_channels')
        return cursor

    def get_precursorquant_headerfields(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT prottable_id '
            'FROM protein_precur_quanted')
        return cursor

    def get_quantchannel_map(self):
        outdict = {}
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT channel_id, prottable_id, channel_name, amount_psms_name'
            ' FROM protquant_channels')
        for channel_id, fnid, channel_name, amount_psms_name in cursor:
            try:
                outdict[fnid][channel_name] = (channel_id, amount_psms_name)
            except KeyError:
                outdict[fnid] = {channel_name: (channel_id, amount_psms_name)}
        return outdict

    def store_isobaric_quants(self, quants):
        self.store_many(
            'INSERT INTO protein_iso_quanted(pacc_id, channel_id, '
            'quantvalue, amount_psms) '
            'VALUES (?, ?, ?, ?)', quants)

    def store_protprob(self, probabilities):
        self.store_many(
            'INSERT INTO protein_probability(pacc_id, prottable_id, '
            'probability) '
            'VALUES (?, ?, ?)', probabilities)
