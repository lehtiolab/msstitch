from app.actions.prottable import precursorarea
from app.dataformats import prottable as prottabledata


def create_protein_table_with_precursor_quant(psms, headerfields):
    top_ms1_psms = precursorarea.generate_top_psms(psms)
    for protacc in top_ms1_psms:
        amount = precursorarea.calculate_protein_precursor_quant(top_ms1_psms,
                                                                 protacc)
        outprotein = {prottabledata.HEADER_PROTEIN: protacc,
                      headerfields['precursorquant'][
                          prottabledata.HEADER_AREA][None]: str(amount)}
        yield outprotein
