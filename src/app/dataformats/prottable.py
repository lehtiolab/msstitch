HEADER_PROTEIN = 'Protein ID'
HEADER_PROTEINS = 'Protein ID(s)'
HEADER_GENEID = 'Gene ID'
HEADER_GENENAME = 'Gene Name'
HEADER_DESCRIPTION = 'Description'
HEADER_COVERAGE = 'Coverage'
HEADER_NO_PROTEIN = 'Protein count'
HEADER_CONTENTPROT = 'Proteins in group'
HEADER_NO_UNIPEP = 'Unique peptide count'
HEADER_NO_PEPTIDE = 'Peptide count'
HEADER_NO_PSM = 'PSM count'
HEADER_NO_PSMS_SUFFIX = ' - Quanted PSM count'
HEADER_NO_FULLQ_PSMS = 'Fully quanted PSM count'
HEADER_AREA = 'MS1 precursor area'
HEADER_QVAL = 'q-value'
HEADER_PEP = 'PEP'
HEADER_QVAL_MODELED = 'q-value (linear model)'
HEADER_QSCORE = 'Q-score best peptide'

PICKED_HEADER = [HEADER_QSCORE, HEADER_QVAL]
ACCESSIONS = {
        'prottable': HEADER_PROTEIN, 
        'genetable': HEADER_GENEID,
        'associdtable': HEADER_GENENAME,
        }
TPROT_HEADER_ACCS = [HEADER_PROTEIN, HEADER_GENENAME, HEADER_GENEID, HEADER_PROTEIN]
