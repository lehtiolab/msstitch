DECOY_PREFIX = 'decoy_'

HEADER_SPECFILE = '#SpecFile'
HEADER_SCANNR = 'ScanNum'
HEADER_SPECSCANID = 'SpecID'
HEADER_CHARGE = 'Charge'
HEADER_PEPTIDE = 'Peptide'
HEADER_PROTEIN = 'Protein'
HEADER_GENE = 'Gene ID'
HEADER_SYMBOL = 'Gene Name'
HEADER_DESCRIPTION = 'Description'
HEADER_PRECURSOR_MZ = 'Precursor'
HEADER_PRECURSOR_QUANT = 'MS1 area'
HEADER_PRECURSOR_FWHM = 'FWHM'
HEADER_MASTER_PROT = 'Master protein(s)'
HEADER_PG_CONTENT = 'Protein group(s) content'
HEADER_PG_AMOUNT_PROTEIN_HITS = 'Amount of matching proteins in group(s)'
HEADER_PG = [HEADER_MASTER_PROT, HEADER_PG_CONTENT,
             HEADER_PG_AMOUNT_PROTEIN_HITS]
HEADER_SVMSCORE = 'percolator svm-score'
HEADER_MISSED_CLEAVAGE = 'missed_cleavage'
HEADER_PSMQ = 'PSM q-value'
HEADER_PEPTIDE_Q = 'peptide q-value'
HEADER_TARGETDECOY = 'TD'
HEADER_MSGFSCORE = 'MSGFScore'
HEADER_EVALUE = 'EValue'
HEADER_SETNAME = 'Biological set'
HEADER_RETENTION_TIME = 'Retention time(min)'
HEADER_INJECTION_TIME = 'Ion injection time(ms)'
HEADER_ION_MOB = 'Ion mobility(Vs/cm2)'
MOREDATA_HEADER = [HEADER_RETENTION_TIME, HEADER_INJECTION_TIME, HEADER_ION_MOB]
PERCO_HEADER = [HEADER_SVMSCORE, HEADER_PSMQ, HEADER_PEPTIDE_Q, HEADER_TARGETDECOY]
