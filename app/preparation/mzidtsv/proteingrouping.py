from hashlib import md5
from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata
import app.sqlite as lookups
from app.preparation.mzidtsv import confidencefilters as conffilt


# FIXME slow performance: around 6 min for 80000 psms of kalamari set. Means
# it will be around 1h for whole set. Which is not great.
# Maybe we can: hit protein group db for quick check if pg exists
# How to check multiple pgs? wildcards?

def get_header_with_proteingroups(header):
    ix = header.index(mzidtsvdata.HEADER_PROTEIN) + 1
    return header[:ix] + mzidtsvdata.HEADER_PG + header[ix:]


def generate_psms_with_proteingroups(fn, oldheader, pgdb, confkey, conflvl,
                                     lower_is_better, unroll=False):
    build_master_db(fn, oldheader, pgdb, confkey, conflvl, lower_is_better,
                    unroll)
    build_content_db(pgdb) 
    ####
    newheader = [mzidtsvdata.HEADER_MASTER_PROT, mzidtsvdata.HEADER_PG_CONTENT,
                 mzidtsvdata.HEADER_PG_AMOUNT_PROTEIN_HITS]
    for i in range(10):
        yield {x: 'speedtest' for x in oldheader + newheader}
    #for line in tsvreader.generate_tsv_psms(fn, oldheader):
    #    pgcontents = [[master] + x for master, x in zip(pgroups.keys(),
    #                                                    pgroups.values())]
    #    psm = {mzidtsvdata.HEADER_MASTER_PROT: ';'.join(pgroups.keys()),
    #           mzidtsvdata.HEADER_PG_CONTENT: ';'.join(
    #               [','.join([str(y) for y in x]) for x in pgcontents]),
    #           mzidtsvdata.HEADER_PG_AMOUNT_PROTEIN_HITS:
    #           count_protein_group_hits(
    #               lineproteins,
    #               pgcontents),
    #           }
    #    psm.update({head: line[head] for head in oldheader})
    #    yield psm


def build_master_db(fn, oldheader, pgdb, confkey, conflvl, lower_is_better,
                    unroll):
    psm_masters = []
    allmasters = {}
    rownr = 0
    for line in tsvreader.generate_tsv_psms(fn, oldheader):
        if not conffilt.passes_filter(line, conflvl, confkey, lower_is_better):
            rownr += 1
            continue
        if unroll:
            lineproteins = get_all_proteins_from_unrolled_psm(line, pgdb)
        else:
            lineproteins = tsvreader.get_proteins_from_psm(line)
        pepprotmap = pgdb.get_protpepmap_from_proteins(lineproteins)
        masters = get_masters(pepprotmap)
        psm_masters.extend([(rownr, x) for x in masters])
        allmasters.update({x: 1 for x in masters})
        rownr += 1
    pgdb.store_masters(allmasters, psm_masters)
    return allmasters


def build_content_db(pgdb):
    protein_groups = [] 
    for master in pgdb.get_all_masters():
        master = master[0]
       	protein_groups.extend(get_protein_group_content(master, pgdb))
    pgdb.store_protein_group_content(protein_groups)
        

def get_protein_group_content(master, pgdb):
    """For each master protein, we generate the protein group proteins
    complete with sequences, psm_ids and scores.
    
    Returns a list of [protein, master, pep_hits, psm_hits, protein_score],
    which is ready to enter the DB table.
    """
    psms = pgdb.get_peptides_from_protein(master)
    protein_group_plus = pgdb.get_proteins_peptides_from_psms(psms)
    proteins_not_in_group = pgdb.filter_proteins_with_missing_peptides(
        [x[0] for x in protein_group_plus], psms)
    pgmap = {}
    for protein, pepseq, score, psm_id in protein_group_plus:
        if protein in proteins_not_in_group:
            continue
        score = int(score)  # MZIDSCORE is INTEGER
        try:
            pgmap[protein][pepseq].append((psm_id, score))
        except KeyError:
            try:
                pgmap[protein][pepseq] = [(psm_id, score)]
            except KeyError:
                pgmap[protein] = {pepseq: [(psm_id, score)]}
    pgmap = [[protein, master, len(peptides), len([psm for psms in peptides.values()
                                           for psm in psms]),
                       sum([psm[1] for psms in peptides.values() for psm in psms])]
             for protein, peptides in pgmap.items()]
    return pgmap
     

def get_masters(ppgraph):
    """From a protein-peptide graph dictionary (keys proteins,
    values peptides), return master proteins aka those which 
    have no proteins whose peptides are supersets of them.
    If a shared master protein is found, report only the first,
    we will sort later. In this case, the master reported here 
    may be temporary."""
    masters = {}
    for protein, peps in ppgraph.items():
        ismaster = True
        peps = set(peps)
        multimaster = set()
        for subprotein, subpeps in ppgraph.items():
            if protein == subprotein:
                continue
            if peps.issubset(subpeps) and peps.union(subpeps) > peps:
                ismaster = False
                break
            elif peps.issubset(subpeps) and peps.intersection(subpeps) == peps:
                multimaster.update({protein, subprotein})
        if ismaster and not multimaster:
            masters[protein] = 1
        elif ismaster and multimaster:
            masters[sorted(list(multimaster))[0]] = 1
    return masters


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

def get_all_proteins_from_unrolled_psm(psm, pgdb):
    psm_id = tsvreader.get_psm_id_from_line(psm)
    return pgdb.get_proteins_for_peptide([psm_id])


def sort_pgroup_peptides(proteins):
    return sort_amounts(proteins, lookups.PEPTIDE_COUNT_INDEX)


def sort_pgroup_psms(proteins):
    return sort_amounts(proteins, lookups.PSM_COUNT_INDEX)


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


def sort_pgroup_score(proteins):
    return sort_amounts(proteins, 5)


def sort_pgroup_coverage(proteins):
    return [proteins]


def sort_alphabet(proteins):
    return [sorted(proteins, key=lambda x: x[2])]
