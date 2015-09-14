from app.dataformats import prottable as prottabledata


def prepare_qvality_input(features, feattype, get_static_percolator_output=False):
    feat_fields = {'probability': prottabledata.HEADER_PROBABILITY,
                   'qvalue': peptabledata.HEADER_QVAL,
                  }
    for feature in features:
        yield feature[feat_fields[feattype]]
