def simple_val_fetch(feature, sqlmap, headerfields, valkey):
    val = feature[sqlmap[valkey]]
    if val is None:
        return {}
    pool = feature[sqlmap['set_name']]
    hfield = headerfields[pool]
    return {hfield: feature[sqlmap[valkey]]}


def fill_mergefeature(outfeature, iso_fun, ms1_fun, prob_fun, fdr_fun, pep_fun,
                      pdata_fun, feature, sqlfieldmap, headerfields,
                      featuredata_map):
    check_feat = {k: v for k, v in outfeature.items()}
    for fun in [iso_fun, ms1_fun, prob_fun, fdr_fun, pep_fun]:
        outfeature.update(fun(feature, sqlfieldmap, headerfields))
    if outfeature == check_feat:
        return
    outfeature.update(pdata_fun(outfeature, featuredata_map, headerfields))
