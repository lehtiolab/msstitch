from app.dataformats import prottable as prottabledata


def prepare_qvality_input(proteins, feattype=False,
                          get_static_percolator_output=False):
    for protein in proteins:
        yield protein[prottabledata.HEADER_PROBABILITY]
