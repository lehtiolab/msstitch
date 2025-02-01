import os
import sys
import re
from app.readers import mzidplus as readers
from app.readers import tsv as tsvreader
from app.readers import xml
from app.readers import percolator as percoreader 
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
            yield scan, specidi
