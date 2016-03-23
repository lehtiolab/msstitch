from app.dataformats import peptable as peptabledata
from app.actions.shared.pepprot_isoquant import base_add_isoquant_data


def add_isoquant_data(peptides, quantpeptides, quantacc, quantfields):
    """Runs through a peptide table and adds quant data from ANOTHER peptide 
    table that contains that data."""
    for peptide in base_add_isoquant_data(peptides, quantpeptides, 
                                          peptabledata.HEADER_PEPTIDE,
                                          quantacc, quantfields):
        yield peptide 
