from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata


def generate_psms_with_proteingroups(fn, oldheader, pgdb, unroll=False):
    for line in tsvreader.generate_tsv_psms_line(fn):
        if unroll:
            lineproteins = tsvreader.get_unrolled_protein(line)
        else:
            lineproteins = tsvreader.get_multiple_proteins(line)
        pgroups = group_proteins(lineproteins, pgdb, unroll)
        pgcontents = pgroups.values()
        psm = {mzidtsvdata.HEADER_MASTER_PROT: ';'.join(pgroups.keys()),
               mzidtsvdata.HEADER_PG_CONTENT: ';'.join(
                   [','.join([str(y) for y in x]) for x in pgcontents]),
               mzidtsvdata.HEADER_PG_AMOUNT_PROTEIN_HITS:
               count_protein_group_hits(
                   lineproteins,
                   pgcontents),
               }
        psm.update({head: value for head, value in zip(oldheader, line)})
        yield psm


def count_protein_group_hits(proteins, pgcontents):
    proteins = set(proteins)
    hit_counts = []
    for pgroup in pgcontents:
        hit_counts.append(len(proteins.intersection(pgroup)))
    return ';'.join(hit_counts)


def group_proteins(proteins, pgdb, unroll=False):
    """Generates protein groups per PSM."""
    pp_graph = get_protpep_graph(proteins, pgdb)
    protein_groups = {x: False for x in pp_graph}
    for protein in pp_graph:
        protein_groups[protein] = get_slave_proteins(protein, pp_graph)
    return {k: v for k, v in protein_groups.items() if v}


def get_protpep_graph(proteins, pgdb):
    """Returns graph as a dict:
        {protein: set([peptide1, peptide2]),}
        This methods calls the db 3 times to get protein groups of the
        proteins passed.
        # TODO See if there is a faster implementation.
    """
    peptides = pgdb.get_peptides_from_proteins(proteins)
    proteingraph = pgdb.get_proteins_peptides_from_peptides(peptides)
    proteins_not_in_group = pgdb.filter_proteins_with_missing_peptides(
        proteingraph.keys(), peptides)
    return {k: v for k, v in proteingraph.items()
            if k not in proteins_not_in_group}


def get_slave_proteins(protein, graph):
    # FIXME ties? other problems?
    slave_proteins = []
    for subprotein, peps in graph.items():
        if subprotein == protein:
            continue
        elif graph[protein].issubset(peps):
            return []
        elif peps.issubset(graph[protein]):
            slave_proteins.append(subprotein)
    return slave_proteins
