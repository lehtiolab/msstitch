from app.lookups.sqlite import proteingroups as lookups


def get_sortfnxs(evidence):
    sortfnxs = [sort_pgroup_peptides,
                sort_pgroup_psms,
                sort_pgroup_score,
                ]
    if evidence:
        sortfnxs.append(sort_evidence_score)
    sortfnxs.extend([sort_pgroup_coverage, sort_alphabet])
    return sortfnxs


def sort_to_get_master(pgroup, evidence):
    sortfnxs = get_sortfnxs(evidence)
    sorted_pg = sort_protein_group(pgroup, sortfnxs, 0)
    return {'master_id': sorted_pg[0][lookups.MASTER_INDEX],
            'protein_acc': sorted_pg[0][lookups.PROTEIN_ACC_INDEX]}


def sort_protein_groups(pgroups, evidence):
    """Gets a protein groups containing dict pgroups, for each master (key)
    there is a list of protein tuples that comprise a proteingroup. This loops
    the groups and returns a sorted group. Assumes master is already sorted
    as the proper master"""
    sortfnxs = get_sortfnxs(evidence)
    pgroups_out = {}
    for master, pgroup in pgroups.items():
        sorted_pg = sort_protein_group(pgroup, sortfnxs, 0)
        pgroups_out[master] = sorted_pg
    return pgroups_out


def sort_protein_group(pgroup, sortfunctions, sortfunc_index):
    """Recursive function that sorts protein group by a number of sorting
    functions."""
    pgroup_out = []
    subgroups = sortfunctions[sortfunc_index](pgroup)
    sortfunc_index += 1
    for subgroup in subgroups:
        if len(subgroup) > 1 and sortfunc_index < len(sortfunctions):
            pgroup_out.extend(sort_protein_group(subgroup,
                                                 sortfunctions,
                                                 sortfunc_index))
        else:
            pgroup_out.extend(subgroup)
    return pgroup_out


def sort_amounts(proteins, sort_index):
    """Generic function for sorting peptides and psms. Assumes a higher
    number is better for what is passed at sort_index position in protein."""
    amounts = {}
    for protein in proteins:
        amount_x_for_protein = protein[sort_index]
        try:
            amounts[amount_x_for_protein].append(protein)
        except KeyError:
            amounts[amount_x_for_protein] = [protein]
    return [v for k, v in sorted(amounts.items(), reverse=True)]


def sort_pgroup_peptides(proteins):
    return sort_amounts(proteins, lookups.PEPTIDE_COUNT_INDEX)


def sort_pgroup_psms(proteins):
    return sort_amounts(proteins, lookups.PSM_COUNT_INDEX)


def sort_pgroup_score(proteins):
    return sort_amounts(proteins, lookups.PROTEIN_SCORE_INDEX)


def sort_pgroup_coverage(proteins):
    return sort_amounts(proteins, lookups.COVERAGE_INDEX)


def sort_evidence_score(proteins):
    return sort_amounts(proteins, lookups.EVIDENCE_LVL_INDEX)


def sort_alphabet(proteins):
    return [sorted(proteins, key=lambda x: x[2])]
