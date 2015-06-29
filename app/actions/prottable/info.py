from app.dataformats import prottable as prottabledata


def add_record_to_proteindata(proteindata, p_acc, pool, psmdata):
    seq, psm_id, desc, cov = psmdata[2], psmdata[3], psmdata[4], psmdata[5]
    try:
        proteindata[p_acc][pool]['psms'].add(psm_id)
    except KeyError:
        emptyinfo = {'psms': set(), 'peptides': set(), 'proteins': set(),
                     'unipeps': 0}
        try:
            proteindata[p_acc][pool] = emptyinfo
        except KeyError:
            proteindata[p_acc] = {pool: emptyinfo, 'desc': desc, 'cov': cov}
    proteindata[p_acc][pool]['psms'].add(psm_id)
    proteindata[p_acc][pool]['peptides'].add(seq)
    #proteindata[p_acc][pool]['proteins'].add(pg_content)


def count_peps_psms(proteindata, p_acc, pool):
    data = proteindata[p_acc][pool]
    proteindata[p_acc][pool]['psms'] = len(data['psms'])
    #proteindata[prot][pool]['proteins'] = len(data['proteins'])
    proteindata[p_acc][pool]['peptides'] = len(data['peptides'])


def create_proteindata_map(pgdb, pool_to_output=False):
    protein_psms_data = pgdb.get_all_protein_psms_with_sets()
    proteindata = {}
    psmdata = next(protein_psms_data)
    last_prot, last_pool = psmdata[0], psmdata[1]
    add_record_to_proteindata(proteindata, last_prot, last_pool, psmdata)
    for psmdata in protein_psms_data:
        p_acc, samplepool = psmdata[0], psmdata[1]
        if pool_to_output and samplepool != pool_to_output:
            continue
        if samplepool != last_pool or p_acc != last_prot:
            count_peps_psms(proteindata, last_prot, last_pool)
            last_pool, last_prot = samplepool, p_acc
        add_record_to_proteindata(proteindata, p_acc, samplepool, psmdata)
    count_peps_psms(proteindata, last_prot, last_pool)
    get_unique_peptides(proteindata)
    return proteindata


def add_protein_data(proteins, pgdb, headerfields, pool_to_output=False):
    """First creates a map with all master proteins with data,
    then outputs protein data dicts for rows of a tsv. If a pool
    is given then only output for that pool will be shown in the
    protein table."""
    proteindata = create_proteindata_map(pgdb, pool_to_output)
    for protein in proteins:
        outprotein = {k: v for k, v in protein.items()}
        protein_acc = protein[prottabledata.HEADER_PROTEIN]
        if not protein_acc in proteindata:
            continue
        outprotein.update(get_protein_data(proteindata, protein_acc,
                                           headerfields))
        outprotein = {k: str(v) for k, v in outprotein.items()}
        yield outprotein


def get_protein_data(proteindata, p_acc, headerfields):
    """Parses protein data for a certain protein into tsv output
    dictionary"""
    proteincount = 'na'
    outdict = {}
    hfields = [prottabledata.HEADER_NO_UNIPEP,
               prottabledata.HEADER_NO_PEPTIDE,
               prottabledata.HEADER_NO_PSM,
               ]
    for pool, pdata in proteindata[p_acc].items():
        pool_values = [pdata['unipeps'], pdata['peptides'], pdata['psms']]
        outdict.update({headerfields['proteindata'][hfield][pool]: val
                        for (hfield, val) in zip(hfields, pool_values)})
    outdict.update({prottabledata.HEADER_DESCRIPTION: pdata[p_acc]['desc'],
                    prottabledata.HEADER_COVERAGE: pdata[p_acc]['cov'],
                    prottabledata.HEADER_NO_PROTEIN: proteincount,
                    })
    return outdict


def get_pep_prot_map(proteindata):
    seq_protein_map = {}
    for protein, pools in proteindata.items():
        for pool, pooldata in pools.items():
            seqs = pooldata['peptides']
        for seq in seqs:
            try:
                seq_protein_map[seq][pool].add(protein)
            except KeyError:
                seq_protein_map[seq][pool] = {protein}
            except KeyError:
                seq_protein_map[seq] = {pool: {protein}}


def get_unique_peptides(proteindata):
    seq_protein_map = get_pep_prot_map(proteindata)
    for seq, pools in seq_protein_map.items():
        for pool, proteins in pools.items():
            if len(proteins) == 1:
                proteindata[next(iter(proteins))][pool]['unipeps'] += 1
