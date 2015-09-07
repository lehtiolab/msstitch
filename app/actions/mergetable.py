def simple_val_fetch(feature, sqlmap, headerfields, poolkey, valkey):
    pool = feature[sqlmap[poolkey]]
    hfield = headerfields[pool]
    return {hfield: feature[sqlmap[valkey]]}


def fill_mergefeature(outfeature, iso_fun, ms1_fun, prob_fun, fdr_fun, pep_fun,
                      pdata_fun, feature, sqlfieldmap, headerfields,
                      featuredata_map):
    for fun in [iso_fun, ms1_fun, prob_fun, fdr_fun, pep_fun]:
        outfeature.update(fun(feature, sqlfieldmap, headerfields))
    outfeature.update(pdata_fun(outfeature, featuredata_map, headerfields))


def parse_NA(feature, header):
    for field in header:
        try:
            feature[field] = str(feature[field])
        except KeyError:
            feature[field] = 'NA'
    return feature
