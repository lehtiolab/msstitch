from app.lookups.sqlite.protpeptable import ProtPepTable


class PepTableProteinCentricDB(ProtPepTable):
    datatype = 'peptide'
    colmap = {'peptide_sequences': ['pep_id', 'sequence'],
              'peptide_precur_quanted': ['pep_id', 'peptable_id', 'quant'],
              'peptide_fdr': ['pep_id', 'peptable_id', 'fdr'],
              'peptide_pep': ['pep_id', 'peptable_id', 'pep'],
              'pepquant_channels': ['channel_id', 'peptable_id',
                                    'channel_name', 'amount_psms_name'],
              'peptide_iso_quanted': ['peptidequant_id', 'pep_id',
                                      'channel_id', 'quantvalue',
                                      'amount_psms'],
              }

    def add_tables(self):
        self.create_tables(['peptide_tables', 'pepquant_channels',
                            'peptide_iso_quanted', 'peptide_precur_quanted',
                            'peptide_fdr', 'peptide_pep'])


class PepTableGeneCentricDB(PepTableProteinCentricDB):
    datatype = 'peptide'

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
        genetable = 'genes'
        return self.get_unique_gene_psms(genetable, fields, firstjoin,
                                         extrajoins)


class PepTablePlainDB(PepTableProteinCentricDB):
    datatype = 'peptide'

    def get_proteins_psms_for_map(self):
        fields = ['p.protein_acc', 'sets.set_name', 'pep.sequence',
                  'psm.psm_id']
        firstjoin = ('protein_psm', 'pp', 'protein_acc')
        return self.get_proteins_psms('proteins', fields, firstjoin)
