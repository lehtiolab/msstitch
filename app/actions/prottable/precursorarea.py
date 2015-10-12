from app.dataformats import prottable as prottabledata
from app.dataformats import mzidtsv as mzidtsvdata


def generate_top_psms(psms, protcol):
    """Fed with a psms generator, this returns the 3 PSMs with
    the highest precursor intensities (or areas, or whatever is
    given in the HEADER_PRECURSOR_QUANT"""
    top_ms1_psms = {}
    for psm in psms:
        protacc = psm[protcol]
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
    return top_ms1_psms


def calculate_protein_precursor_quant(top_ms1_psms, prot_acc):
    try:
        amounts = top_ms1_psms[prot_acc]
    except KeyError:
        return 'NA'
    else:
        amounts = [x for x in amounts if x > 0]
        return sum(amounts) / len(amounts)


def add_ms1_quant_from_top3_mzidtsv(proteins, psms, headerfields, protcol):
    """Collects PSMs with the highes precursor quant values,
    adds sum of the top 3 of these to a protein table"""
    if not protcol:
        protcol = mzidtsvdata.HEADER_MASTER_PROT
    top_ms1_psms = generate_top_psms(psms, protcol)
    for protein in proteins:
        prot_acc = protein[prottabledata.HEADER_PROTEIN]
        prec_area = calculate_protein_precursor_quant(top_ms1_psms, prot_acc)
        outprotein = {k: v for k, v in protein.items()}
        outprotein[headerfields['precursorquant'][
            prottabledata.HEADER_AREA][None]] = str(prec_area)
        yield outprotein
