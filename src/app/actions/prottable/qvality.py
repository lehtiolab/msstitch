from app.dataformats import prottable as prottabledata
from app.dataformats import peptable as peptabledata


def prepare_qvality_input(features, feattype,
                          get_static_percolator_output=False):
    feat_fields = {'probability': prottabledata.HEADER_PROBABILITY,
                   'qvalue': prottabledata.HEADER_QSCORE,
                   'svm': peptabledata.HEADER_SVM_SCORE,
                   }
    for feature in features:
        yield feature[feat_fields[feattype]]
