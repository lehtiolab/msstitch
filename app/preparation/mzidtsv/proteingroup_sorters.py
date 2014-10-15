import app.sqlite as lookups


def sort_protein_groups(pgroups, fasta, evidence_levels):
    sortfnxs = [sort_pgroup_peptides,
                sort_pgroup_psms,
                sort_pgroup_score,
                ]
    if evidence_levels:
        sortfnxs.append(sort_evidence_score)
    if fasta:
        sortfnxs.append(sort_pgroup_coverage)
    sortfnxs.append(sort_alphabet)
    pgroups_out = {}
    for pgroup in pgroups.values():
        sorted_pgroup = sort_protein_group(pgroup, sortfnxs, 0)
        pgroups_out[sorted_pgroup[0][lookups.MASTER_INDEX]] = sorted_pgroup
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
    """Generic function for sorting peptides and psms"""
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
    return [proteins]


def sort_evidence_score(proteins):
    return sort_amounts(proteins, lookups.EVIDENCE_LVL_INDEX)


def sort_coin_toss(proteins):
    return [proteins]


def sort_alphabet(proteins):
    return [sorted(proteins, key=lambda x: x[2])]
