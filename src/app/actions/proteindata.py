def add_protgene_to_protdata(proteindata, p_acc, pgenedata, genecentric, pgcontentmap):
    if genecentric:
        desc = pgenedata[1]
        cov, pgcontent = None, None
        if genecentric == 'genes':
            gene = None
            assoc_id, protid = pgenedata[2], pgenedata[3]
        else:
            assoc_id=None
            gene, protid = pgenedata[2], pgenedata[3]
    else:
        desc, cov = pgenedata[2], pgenedata[3]
        gene, assoc_id = pgenedata[4], pgenedata[5]
        pgcontent = pgcontentmap[p_acc]
        protid = None
    proteindata[p_acc]['desc'] = desc
    proteindata[p_acc]['pgcontent'] = pgcontent
    proteindata[p_acc]['cov'] = cov
    proteindata[p_acc]['gene'].add(gene)
    proteindata[p_acc]['aid'].add(assoc_id)
    proteindata[p_acc]['protein_ids'].add(protid)


def add_psms_to_proteindata(proteindata, p_acc, pool, psmdata):
    """Fill function for create_featuredata_map"""
    seq, psm_id = psmdata[2], psmdata[3]
    try:
        proteindata[p_acc]['pools'][pool]['psms'].add(psm_id)
    except KeyError:
        emptyinfo = {'psms': set(), 'peptides': set(), 'unipeps': 0}
        try:
            proteindata[p_acc]['pools'][pool] = emptyinfo
        except KeyError:
            proteindata[p_acc] = {'pools': {pool: emptyinfo},
                                  'gene': set(), 'aid': set(), 'protein_ids': set()}
        proteindata[p_acc]['pools'][pool]['psms'].add(psm_id)
    proteindata[p_acc]['pools'][pool]['peptides'].add(seq)


def create_featuredata_map(pgdb, psm_fill_fun, pgene_fill_fun, genecentric=False, 
        count_fun=None, pool_to_output=False, get_uniques=False, is_peptides=False):
    """Creates dict of protein data containing PSMs, peptides, proteins in
    protein group, unique peptides, description and coverage. Loops through
    PSM/protein matches and uses a passed fill_fun function to actually
    fill the outputted map.
    """
    if genecentric:
        pgcontentmap = None
    else:
        pgcontentmap = get_proteingroup_content(pgdb)
    protein_psms_data = pgdb.get_proteins_psms_for_map()
    proteindata = {}
    psmdata = next(protein_psms_data)
    last_prot, last_pool = psmdata[0], psmdata[1]
    psm_fill_fun(proteindata, last_prot, last_pool, psmdata)
    for psmdata in protein_psms_data:
        p_acc, samplepool = psmdata[0], psmdata[1]
        if pool_to_output and samplepool != pool_to_output:
            continue
        if samplepool != last_pool or p_acc != last_prot:
            if count_fun is not None:
                count_fun(proteindata, last_prot, last_pool)
            last_pool, last_prot = samplepool, p_acc
        psm_fill_fun(proteindata, p_acc, samplepool, psmdata)
    if count_fun is not None:
        count_fun(proteindata, last_prot, last_pool)
    
    for pgenedata in pgdb.get_protein_gene_symbol_for_map():
        if is_peptides:
            # include sequence (index 1 of record) to fill func
            pgene_fill_fun(proteindata, pgenedata[1], pgenedata[0], 
                    pgenedata, genecentric, pgcontentmap)
        else:
            pgene_fill_fun(proteindata, pgenedata[0], 
                    pgenedata, genecentric, pgcontentmap)
    if get_uniques:
        get_unique_peptides(pgdb, proteindata)
    return proteindata


def get_proteingroup_content(pgdb):
    pgmap = {}
    for master, contentprot in pgdb.get_proteingroup_content():
        try:
            pgmap[master].append(contentprot)
        except KeyError:
            pgmap[master] = [contentprot]
    return pgmap


def get_pep_prot_map(pgdb):
    seq_protein_map = {}
    for p_acc, pool, seq in pgdb.get_all_proteins_psms_for_unipeps():
        try:
            seq_protein_map[pool][seq].add(p_acc)
        except KeyError:
            try:
                seq_protein_map[pool][seq] = {p_acc}
            except KeyError:
                seq_protein_map[pool] = {seq: {p_acc}}
    return seq_protein_map


def get_unique_peptides(pgdb, proteindata):
    seq_protein_map = get_pep_prot_map(pgdb)
    for pool, peptides in seq_protein_map.items():
        for proteins in peptides.values():
            if len(proteins) == 1:
                proteindata[proteins.pop()]['pools'][pool]['unipeps'] += 1
