import re
import os


class PepPSMHeader:
        DECOY_PREFIX = 'decoy_'
        HEADER_GENE = 'Gene ID'
        HEADER_SYMBOL = 'Gene Name'
        HEADER_DESCRIPTION = 'Description'
        HEADER_PRECURSOR_QUANT = 'MS1 area'
        HEADER_PREC_PURITY = 'Precursor ion fraction'
        HEADER_PRECURSOR_FWHM = 'FWHM'
        HEADER_MASTER_PROT = 'Master protein(s)'
        HEADER_PG_CONTENT = 'Protein group(s) content'
        HEADER_PG_AMOUNT_PROTEIN_HITS = 'Total number of matching proteins in group(s)'
        HEADER_PG = [HEADER_MASTER_PROT, HEADER_PG_CONTENT, HEADER_PG_AMOUNT_PROTEIN_HITS]
        HEADER_TARGETDECOY = 'TD'
        HEADER_SETNAME = 'Biological set'
        HEADER_INJECTION_TIME = 'Ion injection time(ms)'
        HEADER_SVMSCORE = 'percolator svm-score'


class PSMTableHeader(PepPSMHeader):
    def __init__(self, header):
        self.header = header

        # detect special cases
        if self.HEADER_SVMSCORE in self.header:
            self.post_percolator()

    def post_percolator(self):
        # if perco is used, we get these. Sometimes it isnt and then there may be others, or not
        # detect this when needed!
        self.HEADER_PSMQ = 'PSM q-value'
        self.HEADER_PSM_PEP = 'PSM PEP'
        self.HEADER_PEPTIDE_Q = 'peptide q-value'
        self.HEADER_PEPTIDE_PEP = 'peptide PEP'
        self.PERCO_HEADER = [self.HEADER_SVMSCORE, self.HEADER_PSMQ, self.HEADER_PSM_PEP,
                self.HEADER_PEPTIDE_Q, self.HEADER_PEPTIDE_PEP, self.HEADER_TARGETDECOY]

    def get_psm_id(self, line, specfncol):
        return '{0}_{1}_{2}'.format(line[specfncol],
                                    line[self.HEADER_SPECSCANID],
                                    line[self.HEADER_PEPTIDE])

    def get_proteins_from_psm(self, line):
        """From a line, return list of proteins reported by Mzid2TSV. When unrolled
        lines are given, this returns the single protein from the line."""
        proteins = line[self.HEADER_PROTEIN].split(';')
        outproteins = []
        for protein in proteins:
            prepost_protein = re.sub('\(pre=.*post=.*\)', '', protein).strip()
            outproteins.append(prepost_protein)
        return outproteins

    def get_psm(self, line, unroll, specfncol):
        """Returns from a PSM line peptide sequence,
        and other information about the PSM.
        Return values:
            specfn          -   str
            psm_id 	 	-   str
            specscanid      -   str
            peptideseq      -   str
            score		-   str
        """
        if specfncol is None:
            specfncol = self.HEADER_SPECFILE
        specfn = line[specfncol]
        psm_id = self.get_psm_id(line, specfncol)
        specscanid = line[self.HEADER_SPECSCANID]
        peptideseq = self.get_psm_sequence(line, unroll)
        score = line[self.HEADER_SCORE]
        return specfn, psm_id, specscanid, peptideseq, score

    def get_psm_sequence(self, line, unroll=False):
        peptideseq = line[self.HEADER_PEPTIDE]
        if unroll and '.' in peptideseq:
            peptideseq = peptideseq.split('.')[1]
        return peptideseq

    def get_pepproteins(self, line, specfncol):
        """Returns from a PSM line peptide sequence,
        and other information about the PSM.
        Return values:
            psm_id          -   str
            proteins        -   list of str
        """
        psm_id = self.get_psm_id(line, specfncol)
        proteins = self.get_proteins_from_psm(line)
        return psm_id, proteins


class MSGFPSMTableHeader(PSMTableHeader):
    id_fields_header = ['MSGFScore', 'EValue']

    def __init__(self, header):
        super().__init__(header)

        self.HEADER_SPECFILE = '#SpecFile'
        self.HEADER_SCANNR = 'ScanNum'
        self.HEADER_SPECSCANID = 'SpecID'
        self.HEADER_CHARGE = 'Charge'
        self.HEADER_PEPTIDE = 'Peptide'
        self.HEADER_PROTEIN = 'Protein'
        # needs adding from msstitch
        # make that depend on this type!
        self.HEADER_ION_MOB = 'Ion mobility(Vs/cm2)'
        self.HEADER_MISSED_CLEAVAGE = 'missed_cleavage'
        self.HEADER_RETENTION_TIME = 'Retention time(min)'
        # exclusives
        self.HEADER_SCORE = 'MSGFScore'
        self.HEADER_EVALUE = 'EValue'
        
        # Post setting SE specific fields:
        self.MOREDATA_HEADER = [self.HEADER_RETENTION_TIME, self.HEADER_INJECTION_TIME,
                self.HEADER_ION_MOB]

    def set_header_with_percolator(self):
        ix = self.header.index(self.HEADER_EVALUE) + 1
        self.post_percolator()
        self.header = self.header[:ix] + self.PERCO_HEADER + self.header[ix:]

    def get_scannr(self, psm):
        return psm[self.HEADER_SCANNR]

    def add_fdr_to_mzidtsv(self, psms, specidis, percodata):
        """Takes PSMs from an mzIdentML and its MSGF+ TSV and a corresponding 
        percolator XML. 
        """
        # mzId results and PSM lines can be zipped
        # but we assume that multi-solution scans (e.g. Isoleucine/leucine)
        # have identical scoring, i.e. MSGF is run with -n1 (default)

        for psm, (scan, specidi) in zip(psms, specidis):
            spfile = os.path.splitext(psm[self.HEADER_SPECFILE])[0]
            psmcharge = psm[self.HEADER_CHARGE]
            perco_id = f'{spfile}_{specidi.attrib["id"]}_{scan}_{psmcharge}_{specidi.attrib["rank"]}'
            try:
                percopsm = percodata[perco_id]
            except KeyError:
                continue
            outpsm = {k: v for k,v in psm.items()}
            psmscores = {
                self.HEADER_SVMSCORE: percopsm['svm'], 
                self.HEADER_PSMQ: percopsm['qval'], 
                self.HEADER_PSM_PEP: percopsm['pep'], 
                }
            try:
                pepscores = {
                    self.HEADER_PEPTIDE_Q: percopsm['pepqval'], 
                    self.HEADER_PEPTIDE_PEP: percopsm['peppep'], 
                    self.HEADER_TARGETDECOY: percopsm['td'],
                    }
            except KeyError:
                pepscores = {
                    self.HEADER_PEPTIDE_Q: 1,
                    self.HEADER_PEPTIDE_PEP: 1,
                    self.HEADER_TARGETDECOY: 1,
                    }
            outpsm.update({**psmscores, **pepscores})
            # Remove all decoy protein matches from target proteins, to ensure downstream
            # processing does not trip up on them, e.g. having a decoy master protein.
            if percopsm['td'] != 'decoy':
                outprots = []
                for prot in outpsm[self.HEADER_PROTEIN].split(';'):
                    if not prot.startswith(self.DECOY_PREFIX):
                        outprots.append(prot)
                outpsm[self.HEADER_PROTEIN] = ';'.join(outprots)
            yield outpsm


class SagePSMTableHeader(PSMTableHeader):
    name = 'sage'
    id_fields_header = ['sage_discriminant_score']

    def __init__(self, header):
        super().__init__(header)

        self.HEADER_SPECFILE = 'filename'
        self.HEADER_SPECSCANID = 'scannr'
        self.HEADER_CHARGE = 'charge'
        self.HEADER_PEPTIDE = 'peptide'
        self.HEADER_PROTEIN = 'proteins'
        self.HEADER_RANK = 'rank'
        self.HEADER_PSM_ID = 'psm_id'
        # not used, but we could
        # self.HEADER_PRECURSOR_MZ = 'expmass'

        self.HEADER_SCORE = 'sage_discriminant_score'

        # is output, need to make sure it wont be added:
        self.HEADER_ION_MOB = 'ion_mobility'
        self.HEADER_RETENTION_TIME = 'rt'
        self.HEADER_MISSED_CLEAVAGE = 'missed_cleavages'

        self.HEADER_PSMQ = 'spectrum_q'
        self.HEADER_PEPTIDE_Q = 'peptide_q'
        self.MS2INT = 'ms2_intensity'

        # Post setting SE specific fields:
        self.MOREDATA_HEADER = [self.HEADER_RETENTION_TIME, self.HEADER_INJECTION_TIME,
                self.HEADER_ION_MOB]

    def set_header_with_percolator(self):
        startix = self.header.index(self.HEADER_SCORE) + 1
        end_ix = self.header.index(self.MS2INT)
        self.post_percolator()
        self.header = self.header[:startix] + self.PERCO_HEADER + self.header[end_ix:]

    def get_scannr(self, psm):
        # controllerType=0 controllerNumber=1 scan=244	
        return psm[self.HEADER_SPECSCANID].split('=')[-1]

    def add_fdr_to_mzidtsv(self, psms, specidis, percodata):
        """Takes PSMs from an mzIdentML and its MSGF+ TSV and a corresponding 
        percolator XML. 
        """
        # we assume that multi-solution scans (e.g. Isoleucine/leucine)
        # have identical scoring, i.e. MSGF is run with -n1 (default)

        print(percodata)
        for psm in psms:
            try:
                percopsm = percodata[psm[self.HEADER_PSM_ID]]
            except KeyError:
                print('no')
                print(psm[self.HEADER_PSM_ID])
                continue
            outpsm = {k: v for k,v in psm.items()}
            psmscores = {
                self.HEADER_SVMSCORE: percopsm['svm'], 
                self.HEADER_PSMQ: percopsm['qval'], 
                self.HEADER_PSM_PEP: percopsm['pep'], 
                }
            try:
                pepscores = {
                    self.HEADER_PEPTIDE_Q: percopsm['pepqval'], 
                    self.HEADER_PEPTIDE_PEP: percopsm['peppep'], 
                    self.HEADER_TARGETDECOY: percopsm['td'],
                    }
            except KeyError:
                pepscores = {
                    self.HEADER_PEPTIDE_Q: 1,
                    self.HEADER_PEPTIDE_PEP: 1,
                    self.HEADER_TARGETDECOY: 1,
                    }
            outpsm.update({**psmscores, **pepscores})
            # Remove all decoy protein matches from target proteins, to ensure downstream
            # processing does not trip up on them, e.g. having a decoy master protein.
            if percopsm['td'] != 'decoy':
                outprots = []
                for prot in outpsm[self.HEADER_PROTEIN].split(';'):
                    if not prot.startswith(self.DECOY_PREFIX):
                        outprots.append(prot)
                outpsm[self.HEADER_PROTEIN] = ';'.join(outprots)
            yield outpsm


def get_psmheader(header):
    for ph in [MSGFPSMTableHeader, SagePSMTableHeader]:
        if all(x in header for x in ph.id_fields_header):
            return ph(header)
    raise RuntimeError('Could not identify what kind of PSM table this is, expect fields '
            'specific to either MSGF+ or Sage')
