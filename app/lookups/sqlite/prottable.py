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

    def get_all_proteins_psms_for_unipeps(self, genecentric):
        fields = ['p.protein_acc', 'sets.set_name',
                  'pep.sequence']
        if genecentric:
            firstjoin = ('protein_psm', 'pp', 'protein_acc')
            firsttable = 'proteins'
        else:
            firstjoin = ('psm_protein_groups', 'ppg', 'master_id')
            firsttable = 'protein_group_master'
        return self.get_proteins_psms(firsttable, fields, firstjoin)

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
               'JOIN protein_tables AS pt ON pt.set_id=bs.set_id'.format(', '.join(selects)))
        sql = self.get_sql_joins_mergetable(sql, joins, 'protein')
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


class GeneTableDB(ProtPepTable):
    datatype = 'gene'
    colmap = {'gene_precur_quanted': ['gene_id', 'genetable_id', 'quant'],
              'gene_fdr': ['gene_id', 'genetable_id', 'fdr'],
              'gene_pep': ['gene_id', 'genetable_id', 'pep'],
              'gene_probability': ['gene_id', 'genetable_id',
                                   'probability'],
              }

    def add_tables(self):
        self.create_tables(['gene_tables', 'gene_iso_quanted',
                            'genequant_channels', 'gene_precur_quanted',
                            'gene_probability', 'gene_fdr',
                            'gene_pep'])

    def store_quant_channels(self, quantchannels):
        self.store_many(
            'INSERT INTO genequant_channels(genetable_id, channel_name, '
            'amount_psms_name) VALUES (?, ?, ?)',
            quantchannels)

    def get_all_proteins_psms_for_unipeps(self, genecentric):
        fields = ['p.gene_acc', 'sets.set_name',
                  'pep.sequence']
        if genecentric:
            firstjoin = ('protein_psm', 'pp', 'protein_acc')
            firsttable = 'genes'
        else:
            firstjoin = ('psm_protein_groups', 'ppg', 'master_id')
            firsttable = 'protein_group_master'
        return self.get_proteins_psms(firsttable, fields, firstjoin)

    def get_proteins_psms_for_map(self):
        """Gets gene-PSM combinations from DB and filters out uniques
        on the fly. Filtering is done since PSM are stored per protein,
        not per gene, so there may be a lot of *plicates"""
        fields = ['p.gene_acc', 'sets.set_name',
                  'pep.sequence', 'psm.psm_id',
                  'pd.description', 'aid.assoc_id']
        firstjoin = ('protein_psm', 'pp', 'protein_acc')
        extrajoins = ('LEFT OUTER JOIN prot_desc AS pd USING(protein_acc) '
                      'LEFT OUTER JOIN associated_ids AS aid USING(protein_acc)'
                      )
        lastgene = None
        gpsms_out, gp_ids = [], []
        for gpsm in self.get_proteins_psms('genes', fields, firstjoin,
                                           extrajoins):
            if gpsm[0] != lastgene:
                for outpsm in gpsms_out:
                    yield outpsm
                lastgene = gpsm[0]
                gpsms_out, gp_ids = [], []
            gp_id = gpsm[0] + gpsm[1] + gpsm[3]
            if gp_id not in gp_ids:
                gp_ids.append(gp_id)
                gpsms_out.append(gpsm)
        for outpsm in gpsms_out:
            yield outpsm

    def prepare_mergetable_sql(self, precursor=False, isobaric=False,
                               probability=False, fdr=False, pep=False):
        selects = ['g.gene_acc', 'bs.set_name']
        selectmap, count = self.update_selects({}, ['p_acc', 'set_name'], 0)
        joins = []
        if isobaric:
            selects.extend(['pc.channel_name',
                            'pc.amount_psms_name', 'giq.quantvalue',
                            'giq.amount_psms'])
            joins.extend([('genequant_channels', 'pc', ['gt']),
                          ('gene_iso_quanted', 'giq', ['g', 'pc'], True),
                          ])
            fld = ['channel', 'isoq_psmsfield', 'isoq_val',
                   'isoq_psms']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if precursor:
            selects.extend(['preq.quant'])
            joins.append(('gene_precur_quanted', 'preq', ['g', 'gt'], True))
            fld = ['preq_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if probability:
            selects.extend(['gprob.probability'])
            joins.append(('gene_probability', 'gprob', ['g', 'gt'], True))
            fld = ['prob_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if fdr:
            selects.extend(['gfdr.fdr'])
            joins.append(('gene_fdr', 'gfdr', ['g', 'gt'], True))
            fld = ['fdr_val']
            selectmap, count = self.update_selects(selectmap, fld, count)
        if pep:
            selects.extend(['gpep.pep'])
            joins.append(('gene_pep', 'gpep', ['g', 'gt'], True))
            fld = ['pep_val']
            selectmap, count = self.update_selects(selectmap, fld, count)

        sql = ('SELECT {} FROM genes AS g JOIN biosets AS bs '
               'JOIN gene_tables AS gt ON gt.set_id=bs.set_id'.format(', '.join(selects)))
        sql = self.get_sql_joins_mergetable(sql, joins, 'gene')
        sql = '{0} ORDER BY g.gene_acc'.format(sql)
        return sql, selectmap

    def get_isoquant_amountpsms_channels(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT channel_name, amount_psms_name '
            'FROM genequant_channels')
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
            'SELECT channel_id, genetable_id, channel_name, amount_psms_name'
            ' FROM genequant_channels')
        for channel_id, fnid, channel_name, amount_psms_name in cursor:
            try:
                outdict[fnid][channel_name] = (channel_id, amount_psms_name)
            except KeyError:
                outdict[fnid] = {channel_name: (channel_id, amount_psms_name)}
        return outdict

    def store_isobaric_quants(self, quants):
        self.store_many(
            'INSERT INTO gene_iso_quanted(gene_id, channel_id, '
            'quantvalue, amount_psms) '
            'VALUES (?, ?, ?, ?)', quants)
