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

    def get_proteins_psms(self, extended=False):
        fields = ['pgm.protein_acc', 'sets.set_name',
                  'pep.sequence']
        joins = [('psm_protein_groups', 'ppg', 'master_id'),
                 ('psms', 'psm', 'psm_id'),
                 ('peptide_sequences', 'pep', 'pep_id'),
                 ('mzml', 'sp', 'spectra_id'),
                 ('mzmlfiles', 'mzfn', 'mzmlfile_id'),
                 ('biosets', 'sets', 'set_id'),
                 ]
        if extended:
            fields.extend(['psm.psm_id', 'pd.description', 'pcov.coverage'])
            joins.extend([('prot_desc', 'pd', 'protein_acc'),
                          ('protein_coverage', 'pcov', 'protein_acc'),
                          ])
        sql = ('SELECT {} FROM protein_group_master '
               'AS pgm'.format(', '.join(fields)))
        join_sql = ' '.join(['JOIN {} AS {} USING({})'.format(
            j[0], j[1], j[2]) for j in joins])
        sql = '{} {} ORDER BY pgm.protein_acc, sets.set_name'.format(sql,
                                                                     join_sql)
        cursor = self.get_cursor()
        return cursor.execute(sql)

    def prepare_mergetable_sql(self, precursor=False, isobaric=False,
                               probability=False, fdr=False, pep=False):
        selects = ['p.protein_acc']
        selectmap, count = self.update_selects({}, ['p_acc'], 0)
        joins = []
        if isobaric:
            selects.extend(['pc.channel_name', 'bs.set_name',
                            'pc.amount_psms_name', 'piq.quantvalue',
                            'piq.amount_psms'])
            joins.extend([('protein_iso_quanted', 'piq', 'p', 'pacc_id', True),
                          ('protquant_channels', 'pc', 'piq', 'channel_id'),
                          ('protein_tables', 'pt', 'pc', 'prottable_id'),
                          ('biosets', 'bs', 'pt', 'set_id')
                          ])
            fld = ['channel', 'isoq_poolname', 'isoq_psmsfield', 'isoq_val',
                   'isoq_psms']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if precursor:
            selects.extend(['prqbs.set_name', 'preq.quant'])
            joins.extend([('protein_precur_quanted', 'preq', 'p', 'pacc_id',
                           True),
                          ('protein_tables', 'prqpt', 'preq', 'prottable_id'),
                          ('biosets', 'prqbs', 'prqpt', 'set_id')
                          ])
            fld = ['preq_poolname', 'preq_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if probability:
            selects.extend(['probbs.set_name', 'pprob.probability'])
            joins.extend([('protein_probability', 'pprob', 'p', 'pacc_id',
                           True),
                          ('protein_tables', 'probpt', 'pprob',
                           'prottable_id'),
                          ('biosets', 'probbs', 'probpt', 'set_id')
                          ])
            fld = ['prob_poolname', 'prob_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if fdr:
            selects.extend(['fdrbs.set_name', 'pfdr.fdr'])
            joins.extend([('protein_fdr', 'pfdr', 'p', 'pacc_id', True),
                          ('protein_tables', 'fdrpt', 'pfdr', 'prottable_id'),
                          ('biosets', 'fdrbs', 'fdrpt', 'set_id')
                          ])
            fld = ['fdr_poolname', 'fdr_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if pep:
            selects.extend(['pepbs.set_name', 'ppep.pep'])
            joins.extend([('protein_pep', 'ppep', 'p', 'pacc_id', True),
                          ('protein_tables', 'peppt', 'ppep', 'prottable_id'),
                          ('biosets', 'pepbs', 'peppt', 'set_id')
                          ])
            fld = ['pep_poolname', 'pep_val']
            selectmap, count = self.update_selects(selectmap, fld, count)

        sql = 'SELECT {} FROM proteins AS p'.format(
            ', '.join(selects))
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
