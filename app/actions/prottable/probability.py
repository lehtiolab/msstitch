from app.dataformats import mzidtsv as mzidtsvdata
from app.dataformats import prottable as prottabledata 


def add_nesvi_protein_probability(proteins, peptides, headerfields):
    protein_probs = {}
    for peptide in peptides:
        protacc = peptide[mzidtsvdata.HEADER_MASTER_PROT]
        pep = peptide[mzidtsvdata.HEADER_PEPTIDE_PEP]
        if ';' in protacc or pep in ['NA', False]:
            continue
        pep = float(pep)
        try:
            protein_probs[protacc] = protein_probs[protacc] * pep
        except KeyError:
            protein_probs[protacc] = pep
    for protein in proteins: 
        outprotein = {k: v for k, v in protein.items()}
        try:
            protprob = protein_probs[protein[prottabledata.HEADER_PROTEIN]]
        except KeyError:
            protprob = 'NA'
        outprotein[headerfields['probability'][prottabledata.HEADER_PROBABILITY][None]] = str(protprob)
        yield outprotein
