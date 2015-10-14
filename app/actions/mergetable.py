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


def parse_NA(feature, header):
    for field in header:
        try:
            feature[field] = str(feature[field])
        except KeyError:
            feature[field] = 'NA'
    return feature


def create_featuredata_map(pgdb, fill_fun, genecentric=False, count_fun=None,
                           pool_to_output=False, get_uniques=False):
    """Creates dict of protein data containing PSMs, peptides, proteins in
    protein group, unique peptides, description and coverage. Loops through
    PSM/protein matches and uses a passed fill_fun function to actually
    fill the outputted map.
    """
    if genecentric:
        pgcontentmap = None
        protein_psms_data = pgdb.get_proteins_psms_genecentric()
    else:
        pgcontentmap = get_proteingroup_content(pgdb)
        protein_psms_data = pgdb.get_proteins_psms_pgrouped()
    proteindata = {}
    psmdata = next(protein_psms_data)
    last_prot, last_pool = psmdata[0], psmdata[1]
    fill_fun(proteindata, last_prot, last_pool, psmdata, genecentric,
             pgcontentmap)
    for psmdata in protein_psms_data:
        p_acc, samplepool = psmdata[0], psmdata[1]
        if pool_to_output and samplepool != pool_to_output:
            continue
        if samplepool != last_pool or p_acc != last_prot:
            if count_fun is not None:
                count_fun(proteindata, last_prot, last_pool)
            last_pool, last_prot = samplepool, p_acc
        fill_fun(proteindata, p_acc, samplepool, psmdata, genecentric,
                 pgcontentmap)
    if count_fun is not None:
        count_fun(proteindata, last_prot, last_pool)
    if get_uniques:
        get_unique_peptides(pgdb, proteindata, genecentric)
    return proteindata


def get_proteingroup_content(pgdb):
    pgmap = {}
    for master, contentprot in pgdb.get_proteingroup_content():
        try:
            pgmap[master].append(contentprot)
        except KeyError:
            pgmap[master] = [contentprot]
    return pgmap


def get_pep_prot_map(pgdb, genecentric):
    seq_protein_map = {}
    for p_acc, pool, seq in pgdb.get_all_proteins_psms_for_unipeps(genecentric):
        try:
            seq_protein_map[pool][seq].add(p_acc)
        except KeyError:
            try:
                seq_protein_map[pool][seq] = {p_acc}
            except KeyError:
                seq_protein_map[pool] = {seq: {p_acc}}
    return seq_protein_map


def get_unique_peptides(pgdb, proteindata, genecentric):
    seq_protein_map = get_pep_prot_map(pgdb, genecentric)
    for pool, peptides in seq_protein_map.items():
        for proteins in peptides.values():
            if len(proteins) == 1:
                proteindata[
                    next(iter(proteins))]['pools'][pool]['unipeps'] += 1
