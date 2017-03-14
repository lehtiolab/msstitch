from app.readers import tsv as reader


def base_add_isoquant_data(features, quantfeatures, acc_col, quantacc_col,
                           quantfields):
    """Generic function that takes a peptide or protein table and adds
    quant data from ANOTHER such table."""
    quant_map = get_quantmap(quantfeatures, quantacc_col, quantfields)
    for feature in features:
        feat_acc = feature[acc_col]
        outfeat = {k: v for k, v in feature.items()}
        try:
            outfeat.update(quant_map[feat_acc])
        except KeyError:
            outfeat.update({field: 'NA' for field in quantfields})
        yield outfeat


def get_quantmap(features, acc_col, quantfields):
    """Runs through proteins that are in a quanted protein table, extracts
    and maps their information based on the quantfields list input.
    Map is a dict with protein_accessions as keys."""
    qmap = {}
    for feature in features:
        feat_acc = feature.pop(acc_col)
        qmap[feat_acc] = {qf: feature[qf] for qf in quantfields}
    return qmap

