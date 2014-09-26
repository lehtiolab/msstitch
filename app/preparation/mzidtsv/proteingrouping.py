from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata
from app.sqlite import ProteinPeptideDB
from app.preparation.mzidtsv import confidencefilters as conffilt


def get_header_with_proteingroups(header):
    ix = header.index(mzidtsvdata.HEADER_PROTEIN) + 1
    return header[:ix] + mzidtsvdata.HEADER_PG + header[ix:]


def generate_psms_with_proteingroups(fn, oldheader, pgdbfn, confkey, conflvl,
                                     lower_is_better, unroll=False):
    pgdb = ProteinPeptideDB(pgdbfn)
    for line in tsvreader.generate_tsv_psms(fn, oldheader):
        if not conffilt.passes_filter(line, conflvl, confkey, lower_is_better):
            continue
        if unroll:
            lineproteins = get_all_proteins_from_unrolled_psm(line, pgdb)
        else:
            lineproteins = tsvreader.get_proteins_from_psm(line)
        pgroups = group_proteins(lineproteins, pgdb)
        pgcontents = [[master] + x for master, x in zip(pgroups.keys(),
                                                        pgroups.values())]
        psm = {mzidtsvdata.HEADER_MASTER_PROT: ';'.join(pgroups.keys()),
               mzidtsvdata.HEADER_PG_CONTENT: ';'.join(
                   [','.join([str(y) for y in x]) for x in pgcontents]),
               mzidtsvdata.HEADER_PG_AMOUNT_PROTEIN_HITS:
               count_protein_group_hits(
                   lineproteins,
                   pgcontents),
               }
        psm.update({head: line[head] for head in oldheader})
        yield psm


def count_protein_group_hits(proteins, pgcontents):
    proteins = set(proteins)
    hit_counts = []
    for pgroup in pgcontents:
        hit_counts.append(str(len(proteins.intersection(pgroup))))
    return ';'.join(hit_counts)


def group_proteins(proteins, pgdb):
    """Generates protein groups per PSM. First calls get_slave_proteins to
    determine which of a graph is the master protein. Then calls a sorting
    function to sort the rest of the group."""
    sortfnxs = [sort_pgroup_peptides,
                sort_pgroup_psms,
                sort_pgroup_score,
                sort_pgroup_coverage,
                sort_alphabet,
                ]
    pp_graph = get_protpep_graph(proteins, pgdb)
    protein_groups = {}
    for protein in pp_graph:
        slaves = get_slave_proteins(protein, pp_graph)
        if slaves is not False:
            protein_groups[protein] = slaves
    for protein, pgroup in protein_groups.items():
        protein_groups[protein] = sort_proteingroup(sortfnxs, 0,
                                                    pgroup, pp_graph)
    return protein_groups


def sort_proteingroup(sortfunctions, sortfunc_index, pgroup, ppgraph):
    """Recursive function that sorts protein group by a number of sorting
    functions."""
    pgroup_out = []
    subgroups = sortfunctions[sortfunc_index](pgroup, ppgraph)
    sortfunc_index += 1
    for subgroup in subgroups:
        if len(subgroup) > 1 and sortfunc_index < len(sortfunctions):
            pgroup_out.extend(sort_proteingroup(sortfunctions,
                                                sortfunc_index,
                                                subgroup, ppgraph))
        else:
            pgroup_out.extend(subgroup)
    return pgroup_out


def get_all_proteins_from_unrolled_psm(psm, pgdb):
    psm_id = tsvreader.get_psm_id_from_line(psm)
    return pgdb.get_proteins_for_peptide([psm_id])


def get_protpep_graph(proteins, pgdb):
    """Returns graph as a dict:
        {protein: {peptide1: [(psm_id, score), (psm_id, score)],}}
        This methods calls the db 3 times to get protein groups of the
        proteins passed.
        # TODO See if there is a faster implementation.
    """
    peptides = pgdb.get_peptides_from_proteins(proteins)
    proteingraph = pgdb.get_proteins_peptides_from_peptides(peptides)
    proteins_not_in_group = pgdb.filter_proteins_with_missing_peptides(
        list(proteingraph.keys()), peptides)
    return {k: v for k, v in proteingraph.items()
            if k not in proteins_not_in_group}


def sort_pgroup_peptides(proteins, ppgraph):
    return sort_evaluator(sort_amounts(proteins, ppgraph))


def sort_pgroup_psms(proteins, ppgraph):
    return sort_evaluator(sort_amounts(proteins, ppgraph, innerlookup=True))


def sort_amounts(proteins, ppgraph, innerlookup=False):
    """Generic function for sorting peptides and psms"""
    amounts = {}
    for protein in proteins:
        if not innerlookup:
            amount_x_for_protein = len(ppgraph[protein])
        else:
            amount_x_for_protein = len([x for y in ppgraph[protein].values()
                                        for x in y])
        try:
            amounts[amount_x_for_protein].append(protein)
        except KeyError:
            amounts[amount_x_for_protein] = [protein]
    return amounts


def sort_pgroup_score(proteins, ppgraph):
    scores = {}
    for protein in proteins:
        protein_score = sum([int(x[1]) for y in ppgraph[protein].values() for x in y])
        try:
            scores[protein_score].append(protein)
        except KeyError:
            scores[protein_score] = [protein]
    return sort_evaluator(scores)


def sort_evaluator(sort_dict):
    return [v for k, v in sorted(sort_dict.items(), reverse=True)]


def sort_pgroup_coverage(proteins, ppgraph):
    return [proteins]


def sort_alphabet(proteins, ppgraph):
    return [sorted(proteins)]
# FIXME sequence coverage, we need database for that


def get_slave_proteins(protein, graph):
    # FIXME ties? other problems?
    slave_proteins = []
    for subprotein, peps in graph.items():
        if subprotein == protein:
            continue
        elif set(graph[protein]).issubset(peps):
            return False
        elif set(peps).issubset(graph[protein]):
            slave_proteins.append(subprotein)
    return slave_proteins
