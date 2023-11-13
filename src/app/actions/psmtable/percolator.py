import os
import sys
import re
from app.readers import mzidplus as readers
from app.readers import tsv as tsvreader
from app.readers import xml
from app.readers import percolator as percoreader 
from app.dataformats import mzidtsv as psmheaders
from app.dataformats import qvality as qvalityheaders


def get_scores_fdrs_from_percoqvality(percofn, qvality_psms, qvality_peps):
    '''From percolator output get svm scores, and use a separate qvality output
    for getting the calculated FDR values. Nice in case of recalculating FDR
    after filtering etc
    '''
    ns = xml.get_namespace_from_top(percofn, None)
    psms = {}
    for qvalitypsm, psm in zip(qvality_psms, percoreader.generate_psms(percofn, ns)):
        psms[psm.attrib['{%s}psm_id' % ns['xmlns']]] = {
                'svm': float(psm.find('{%s}svm_score' % ns['xmlns']).text),
                'qval': float(qvalitypsm[qvalityheaders.HEADER_QVAL]),
                'pep': float(qvalitypsm[qvalityheaders.HEADER_PEP]),
                'td': 'decoy' if psm.attrib['{%s}decoy' % ns['xmlns']] == 'true' else 'target',
                }
    for qvalitypep, peptide in zip(qvality_peps, percoreader.generate_peptides(percofn, ns)):
        psm_ids = peptide.find('{%s}psm_ids' % ns['xmlns'])
        for psm_id in psm_ids:
            try:
                psms[psm_id.text].update({
                    'pepqval': float(qvalitypep[qvalityheaders.HEADER_QVAL]),
                    'peppep': float(qvalitypep[qvalityheaders.HEADER_PEP]),
                    })
            except KeyError:
                # Edgecase, sometimes filtering percolator data causes this:
                # when there is a peptide with a PSM ID not in percolator PSMs 
                # set all its PSMs pep-q-values to 1
                [psms[x.text].update({'pepqval': 1, 'peppep': 1}) for x in psm_ids if x.text in psms]
                break
    return psms



def get_scores_fdrs_from_percoxml(percofn):
    """
    From a percolator XML output file, take the FDR (q-value) and PEP and apply
    to the PSM table
    """
    ns = xml.get_namespace_from_top(percofn, None)
    # Get PSM FDRs
    psms = {
            p.attrib['{%s}psm_id' % ns['xmlns']]: {
                'svm': float(p.find('{%s}svm_score' % ns['xmlns']).text),
                'qval': float(p.find('{%s}q_value' % ns['xmlns']).text),
                'pep': float(p.find('{%s}pep' % ns['xmlns']).text),
                'td': 'decoy' if p.attrib['{%s}decoy' % ns['xmlns']] == 'true' else 'target',
                }
            for p in percoreader.generate_psms(percofn, ns)
            }
    # Then get peptide FDRs and add them to PSMs
    for peptide in percoreader.generate_peptides(percofn, ns):
        psm_ids = peptide.find('{%s}psm_ids' % ns['xmlns'])
        for psm_id in psm_ids:
            try:
                psms[psm_id.text].update({
                    'pepqval': float(peptide.find('{%s}q_value' % ns['xmlns']).text),
                    'peppep': float(peptide.find('{%s}pep' % ns['xmlns']).text),
                    })
            except KeyError:
                # Edgecase, sometimes filtering percolator data causes this:
                # when there is a peptide with a PSM ID not in percolator PSMs 
                # set all its PSMs pep-q-values to 1
                [psms[x.text].update({'pepqval': 1, 'peppep': 1}) for x in psm_ids if x.text in psms]
                break
    return psms


def generate_spec_id_items(mzid_specidr, mzns):
    scan = 0
    has_scans = True
    for specidr in mzid_specidr:
        specid = specidr.attrib['spectrumID']
        if has_scans:
            try:
                scan = int({x.split('=')[0]: x.split('=')[1] for x in specidr.attrib['spectrumID'].split(' ')}['scan'])
            except KeyError:
                # in e.g. timstof data there are no true scan numbers, percolator sets it by increment
                has_scans = False
                scan += 1
        else:
            scan += 1
        for specidi in specidr.findall('{%s}SpectrumIdentificationItem' % mzns['xmlns']):
            yield scan, specid, specidi


def add_fdr_to_mzidtsv(psms, mzid_specidr, mzns, percodata):
    """Takes PSMs from an mzIdentML and its MSGF+ TSV and a corresponding 
    percolator XML. 
    """
    # mzId results and PSM lines can be zipped
    # but we assume that multi-solution scans (e.g. Isoleucine/leucine)
    # have identical scoring, i.e. MSGF is run with -n1 (default)

    specidis = generate_spec_id_items(mzid_specidr, mzns)
    for psm in psms:
        psmspecid = psm[psmheaders.HEADER_SPECSCANID]
        psmseq = re.sub('[^A-Z]', '', psm[psmheaders.HEADER_PEPTIDE].upper())
        psmcharge = psm[psmheaders.HEADER_CHARGE]

        mzid_specid = False
        spfile = os.path.splitext(psm[psmheaders.HEADER_SPECFILE])[0]
        percopsm = False
        while mzid_specid != psmspecid:
            try:
                scan, mzid_specid, specidi = next(specidis)
            except StopIteration:
                print('PSM tsv file and corresponding mzIdentML are not lined '
                        'up properly')
                sys,exit(1)
            else:
                perco_id = f'{spfile}_{specidi.attrib["id"]}_{scan}_{psmcharge}_{specidi.attrib["rank"]}'
                try:
                    percopsm = percodata[perco_id]
                except KeyError:
                    pass
        if not percopsm:
            continue
        outpsm = {k: v for k,v in psm.items()}
        psmscores = {
            psmheaders.HEADER_SVMSCORE: percopsm['svm'], 
            psmheaders.HEADER_PSMQ: percopsm['qval'], 
            psmheaders.HEADER_PSM_PEP: percopsm['pep'], 
            }
        try:
            pepscores = {
                psmheaders.HEADER_PEPTIDE_Q: percopsm['pepqval'], 
                psmheaders.HEADER_PEPTIDE_PEP: percopsm['peppep'], 
                psmheaders.HEADER_TARGETDECOY: percopsm['td'],
                }
        except KeyError:
            pepscores = {
                psmheaders.HEADER_PEPTIDE_Q: 1,
                psmheaders.HEADER_PEPTIDE_PEP: 1,
                psmheaders.HEADER_TARGETDECOY: 1,
                }
        outpsm.update({**psmscores, **pepscores})
        # Remove all decoy protein matches from target proteins, to ensure downstream
        # processing does not trip up on them, e.g. having a decoy master protein.
        if percopsm['td'] != 'decoy':
            outprots = []
            for prot in outpsm[psmheaders.HEADER_PROTEIN].split(';'):
                if not prot.startswith(psmheaders.DECOY_PREFIX):
                    outprots.append(prot)
            outpsm[psmheaders.HEADER_PROTEIN] = ';'.join(outprots)
        yield outpsm


def get_header_with_percolator(oldheader):
    ix = oldheader.index(psmheaders.HEADER_EVALUE) + 1
    return oldheader[:ix] + psmheaders.PERCO_HEADER + oldheader[ix:]
