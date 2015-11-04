from app.lookups.sqlite.protpeptable import ProtPepTable


class ProtTableDB(ProtPepTable):
    datatype = 'protein'
    colmap = {'proteins': ['pacc_id', 'protein_acc'],
              'protein_precur_quanted': ['pacc_id', 'prottable_id', 'quant'],
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
        # FIXME isnt genecentric ready to removed since the DB interface
        # also changes with genecentric/not gene centric?
        fields = ['p.protein_acc', 'sets.set_name',
                  'pep.sequence']
        if genecentric:
            firstjoin = ('protein_psm', 'pp', 'protein_acc')
            firsttable = 'proteins'
        else:
            firstjoin = ('psm_protein_groups', 'ppg', 'master_id')
            firsttable = 'protein_group_master'
        return self.get_proteins_psms(firsttable, fields, firstjoin)

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
    colmap = {'genes': ['gene_id', 'gene_acc', 'protein_acc'],
              'gene_precur_quanted': ['gene_id', 'genetable_id', 'quant'],
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
        # FIXME isnt genecentric ready to removed since the DB interface
        # also changes with genecentric/not gene centric?
        # I mean: gene table of which the first table is protein_group_master?
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
                      'LEFT OUTER JOIN associated_ids AS aid '
                      'USING(protein_acc)'
                      )
        return self.get_unique_gene_psms(fields, firstjoin, extrajoins)

    def get_unique_gene_psms(self, fields, firstjoin, extrajoins):
        genetable = self.table_map[self.datatype]['feattable']
        lastgene = None
        gpsms_out, gp_ids = [], []
        for gpsm in self.get_proteins_psms(genetable, fields, firstjoin,
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
            'FROM {}'.format(self.table_map[self.datatype]['prectable']))
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
            'INSERT INTO {}(gene_id, channel_id, quantvalue, amount_psms) '
            'VALUES '
            '(?, ?, ?, ?)'.format(self.table_map[self.datatype]['isoqtable']),
            quants)


class GeneTableAssocIDsDB(GeneTableDB):
    datatype = 'assoc'

    def add_tables(self):
        self.colmap.pop('genes')
        self.colmap = {table.replace('gene', 'assoc'): cols
                       for table, cols in self.colmap.items()}
        self.colmap['associated_ids'] = ['gene_id', 'assoc_id', 'protein_acc']
        self.create_tables(['gene_tables', 'assoc_iso_quanted',
                            'genequant_channels', 'assoc_precur_quanted',
                            'assoc_probability', 'assoc_fdr',
                            'assoc_pep'])

    def get_all_proteins_psms_for_unipeps(self, genecentric):
        # FIXME isnt genecentric ready to removed since the DB interface
        # also changes with genecentric/not gene centric?
        # this is a test function and it is not used at all.
        fields = ['p.assoc_id', 'sets.set_name',
                  'pep.sequence']
        firstjoin = ('protein_psm', 'pp', 'protein_acc')
        firsttable = 'associated_ids'
        return self.get_proteins_psms(firsttable, fields, firstjoin)

    def get_proteins_psms_for_map(self):
        """Gets gene-PSM combinations from DB and filters out uniques
        on the fly. Filtering is done since PSM are stored per protein,
        not per gene, so there may be a lot of *plicates"""
        fields = ['p.assoc_id', 'sets.set_name',
                  'pep.sequence', 'psm.psm_id',
                  'pd.description', 'g.gene_acc']
        firstjoin = ('protein_psm', 'pp', 'protein_acc')
        extrajoins = ('LEFT OUTER JOIN prot_desc AS pd USING(protein_acc) '
                      'LEFT OUTER JOIN genes AS g '
                      'USING(protein_acc)'
                      )
        return self.get_unique_gene_psms(fields, firstjoin, extrajoins)
