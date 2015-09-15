from app.dataformats import prottable as prottabledata


def prepare_qvality_input(features, feattype, get_static_percolator_output=False):
    feat_fields = {'probability': prottabledata.HEADER_PROBABILITY,
                   'qvalue': prottabledata.HEADER_BEST_PEPTIDE_Q,
                  }
    for feature in features:
        yield feature[feat_fields[feattype]]
