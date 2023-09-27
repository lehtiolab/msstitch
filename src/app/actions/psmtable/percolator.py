import os
from app.readers import mzidplus as readers
from app.readers import tsv as tsvreader
from app.readers import xml
from app.readers import percolator as percoreader 
from app.dataformats import mzidtsv as psmheaders


def calculate_target_decoy_competition(percofn):
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


def add_fdr_to_mzidtsv(psms, mzid_specidr, mzns, percodata):
    """Takes PSMs from an mzIdentML and its MSGF+ TSV and a corresponding 
    percolator XML. 
    """
    # mzId results and PSM lines can be zipped
    scan = 0
    for specidr in mzid_specidr:
        try:
            scan = int({x.split('=')[0]: x.split('=')[1] for x in specidr.attrib['spectrumID'].split(' ')}['scan'])
        except KeyError:
            # in e.g. timstof data there are no true scan numbers, percolator sets it by increment
            scan += 1
        for specidi in specidr.findall('{%s}SpectrumIdentificationItem' % mzns['xmlns']):
            psm = next(psms)
            # percolator psm ID is: samplename_SII_scanindex_sii.ix_scannr_charge_rank
            # SII_scanindex_sii.ix is in SII.attrib['id']
            spfile = os.path.splitext(psm[psmheaders.HEADER_SPECFILE])[0]
            try:
                percopsm = percodata['{fn}_{sii_id}_{sc}_{ch}_{rk}'.format(fn=spfile,
                    sii_id=specidi.attrib['id'], sc=scan, rk=specidi.attrib['rank'], ch=psm['Charge'])]
            except KeyError:
                continue
            outpsm = {k: v for k,v in psm.items()}
            outpsm.update({
                psmheaders.HEADER_SVMSCORE: percopsm['svm'], 
                psmheaders.HEADER_PSMQ: percopsm['qval'], 
                psmheaders.HEADER_PSM_PEP: percopsm['pep'], 
                psmheaders.HEADER_PEPTIDE_Q: percopsm['pepqval'], 
                psmheaders.HEADER_PEPTIDE_PEP: percopsm['peppep'], 
                psmheaders.HEADER_TARGETDECOY: percopsm['td'],
                })
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
