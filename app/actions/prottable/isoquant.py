from app.dataformats import prottable as prottabledata
from app.actions.shared.pepprot_isoquant import base_add_isoquant_data


def add_isoquant_data(proteins, quantproteins, quantacc, quantfields):
    """Runs through a protein table and adds quant data from ANOTHER protein
    table that contains that data."""
    for protein in base_add_isoquant_data(proteins, quantproteins, 
                                          prottabledata.HEADER_PROTEIN,
                                          quantacc, quantfields):
        yield protein
