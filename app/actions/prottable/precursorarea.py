from app.dataformats import prottable as prottabledata
from app.dataformats import mzidtsv as mzidtsvdata


def add_ms1_quant_from_top3_mzidtsv(proteins, psms, headerfields):
    """Collects PSMs with the highes precursor quant values,
    adds sum of the top 3 of these to a protein table"""
    top_ms1_psms = {}
    for psm in psms:
        protacc = psm[mzidtsvdata.HEADER_MASTER_PROT]
        precursor_amount = psm[mzidtsvdata.HEADER_PRECURSOR_QUANT]
        if ';' in protacc or precursor_amount == 'NA':
            continue
        precursor_amount = float(precursor_amount)
        psm_seq = psm[mzidtsvdata.HEADER_PEPTIDE]
        try:
            min_precursor_amount = min(top_ms1_psms[protacc])
        except KeyError:
            top_ms1_psms[protacc] = {-1: None, -2: None,
                                         precursor_amount: psm_seq}
            continue
        else:
            if precursor_amount > min_precursor_amount:
                top_ms1_psms[protacc][precursor_amount] = psm_seq 
                top_ms1_psms[protacc].pop(min_precursor_amount)
    for protein in proteins:
        outprotein = {k: v for k, v in protein.items()}
        try:
            amounts = top_ms1_psms[protein[prottabledata.HEADER_PROTEIN]]
        except KeyError:
            prec_area = 'NA'
        else:
            amounts = [x for x in amounts if x > 0]
            prec_area = sum(amounts) / len(amounts)
        outprotein[headerfields['precursorquant'][prottabledata.HEADER_AREA][None]] = str(prec_area)
        yield outprotein
