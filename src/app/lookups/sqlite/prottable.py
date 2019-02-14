from app.lookups.sqlite.protpeptable import ProtPepTable


class ProtTableDB(ProtPepTable):
    datatype = 'protein'
    colmap = {'protein_group_master': ['master_id', 'protein_acc'],
              'proteins': ['pacc_id', 'protein_acc'],
              'protein_precur_quanted': ['pacc_id', 'prottable_id', 'quant'],
              'protein_fdr': ['pacc_id', 'prottable_id', 'fdr'],
              'protein_pep': ['pacc_id', 'prottable_id', 'pep'],
              'protein_probability': ['pacc_id', 'prottable_id',
                                      'probability'],
              'protquant_channels': ['channel_id', 'prottable_id',
                                     'channel_name', 'amount_psms_name'],
              'protein_iso_quanted': ['proteinquant_id', 'pacc_id',
                                      'channel_id', 'quantvalue',
                                      'amount_psms'],
              }

    def add_tables(self):
        self.create_tables(['protein_tables', 'protein_iso_quanted',
                            'protquant_channels', 'protein_precur_quanted',
                            'protein_probability', 'protein_fdr',
                            'protein_pep'])

    def get_unique_peptide_nrs(self):
        return self.get_unique_peptide_nrs_base('protein_acc', 'JOIN protein_group_content USING(protein_acc) JOIN protein_group_master AS acctable USING(master_id)')

    def get_precursorquant_headerfields(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT prottable_id '
            'FROM protein_precur_quanted')
        return cursor


class GeneTableDB(ProtPepTable):
    datatype = 'gene'
    colmap = {'genes': ['gene_id', 'gene_acc', 'protein_acc'],
              'gene_precur_quanted': ['gene_id', 'genetable_id', 'quant'],
              'gene_fdr': ['gene_id', 'genetable_id', 'fdr'],
              'gene_pep': ['gene_id', 'genetable_id', 'pep'],
              'gene_probability': ['gene_id', 'genetable_id',
                                   'probability'],
              'genequant_channels': ['channel_id', 'genetable_id',
                                     'channel_name', 'amount_psms_name'],
              'gene_iso_quanted': ['genequant_id', 'gene_id',
                                   'channel_id', 'quantvalue', 'amount_psms'],
              }

    def add_tables(self):
        self.create_tables(['gene_tables', 'gene_iso_quanted',
                            'genequant_channels', 'gene_precur_quanted',
                            'gene_probability', 'gene_fdr',
                            'gene_pep'])

    def get_protein_gene_symbol_for_map(self):
        fields = ['p.gene_acc', 'pd.description',
                  'aid.assoc_id', 'prot.protein_acc']
        genetable = self.table_map[self.datatype]['feattable']
        sql = (
                'SELECT {} FROM {} AS p '
                'LEFT OUTER JOIN prot_desc AS pd USING(protein_acc) '
                'LEFT OUTER JOIN associated_ids AS aid USING(protein_acc) '
                'LEFT OUTER JOIN proteins AS prot USING(protein_acc)')
        cursor = self.get_cursor()
        sql = sql.format(','.join(fields), genetable)
        return cursor.execute(sql)

    def get_proteins_psms_for_map(self):
        """Gets gene-PSM combinations from DB and filters out uniques
        on the fly. Filtering is done since PSM are stored per protein,
        not per gene, so there may be a lot of *plicates"""
        fields = ['p.gene_acc', 'sets.set_name',
                  'pep.sequence', 'psm.psm_id']
        firstjoin = ('protein_psm', 'pp', 'protein_acc')
        genetable = self.table_map[self.datatype]['feattable']
        return self.get_unique_gene_psms(genetable, fields, firstjoin)
             
    def get_unique_peptide_nrs(self):
        return self.get_unique_peptide_nrs_base('gene_acc', 'JOIN genes AS acctable USING(protein_acc)')

    def get_precursorquant_headerfields(self):
        cursor = self.get_cursor()
        cursor.execute(
            'SELECT DISTINCT prottable_id '
            'FROM {}'.format(self.table_map[self.datatype]['prectable']))
        return cursor


class GeneTableAssocIDsDB(GeneTableDB):
    datatype = 'assoc'

    def __init__(self, fn=None):
        super().__init__(fn)
        self.colmap.pop('genes')
        self.colmap = {table.replace('gene', 'assoc'): cols
                       for table, cols in self.colmap.items()}
        self.colmap['genequant_channels'] = self.colmap.pop(
            'assocquant_channels')
        self.colmap['associated_ids'] = ['gene_id', 'assoc_id', 'protein_acc']

    def add_tables(self):
        self.create_tables(['gene_tables', 'assoc_iso_quanted',
                            'genequant_channels', 'assoc_precur_quanted',
                            'assoc_probability', 'assoc_fdr',
                            'assoc_pep'])

    def get_unique_peptide_nrs(self):
        return self.get_unique_peptide_nrs_base('assoc_id', 'JOIN associated_ids AS acctable USING(protein_acc)')

    def get_protein_gene_symbol_for_map(self):
        fields = ['p.assoc_id', 'pd.description',
                  'gene.gene_acc', 'prot.protein_acc']
        genetable = self.table_map[self.datatype]['feattable']
        sql = (
                'SELECT {} FROM {} AS p '
                'LEFT OUTER JOIN prot_desc AS pd USING(protein_acc) '
                'LEFT OUTER JOIN genes AS gene USING(protein_acc) '
                'LEFT OUTER JOIN proteins AS prot USING(protein_acc)')
        cursor = self.get_cursor()
        sql = sql.format(','.join(fields), genetable)
        return cursor.execute(sql)

    def get_proteins_psms_for_map(self):
        """Gets gene-PSM combinations from DB and filters out uniques
        on the fly. Filtering is done since PSM are stored per protein,
        not per gene, so there may be a lot of *plicates"""
        fields = ['p.assoc_id', 'sets.set_name',
                  'pep.sequence', 'psm.psm_id']
        firstjoin = ('protein_psm', 'pp', 'protein_acc')
        genetable = self.table_map[self.datatype]['feattable']
        return self.get_unique_gene_psms(genetable, fields, firstjoin)
