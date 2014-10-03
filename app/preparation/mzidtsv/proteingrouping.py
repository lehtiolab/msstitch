from hashlib import md5
from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata
import app.sqlite as lookups
from app.preparation.mzidtsv import confidencefilters as conffilt


def get_header_with_proteingroups(header):
    ix = header.index(mzidtsvdata.HEADER_PROTEIN) + 1
    return header[:ix] + mzidtsvdata.HEADER_PG + header[ix:]


def get_all_proteins_from_unrolled_psm(psm, pgdb):
    psm_id = tsvreader.get_psm_id_from_line(psm)
    return pgdb.get_proteins_for_peptide([psm_id])


def generate_psms_with_proteingroups(fn, oldheader, newheader, pgdb, confkey, conflvl,
                                     lower_is_better, unroll=False):
    build_master_db(fn, oldheader, pgdb, confkey, conflvl, lower_is_better,
                    unroll)
    build_content_db(pgdb) 
    rownr = 0
    all_protein_group_content = pgdb.get_all_psms_proteingroups()
    protein = next(all_protein_group_content)
    for psm in tsvreader.generate_tsv_psms(fn, oldheader):
        #if rownr % 10000 == 0 and rownr != 0:
        #    break
        if not conffilt.passes_filter(psm, conflvl, confkey, lower_is_better):
            rownr += 1
            continue
        if unroll:
            lineproteins = get_all_proteins_from_unrolled_psm(psm, pgdb)
        else:
            lineproteins = tsvreader.get_proteins_from_psm(psm)
        proteins_in_groups = {}
        while protein[0] == rownr:
            try:
                proteins_in_groups[protein[lookups.MASTER_INDEX]].append(protein)
            except KeyError:
                proteins_in_groups[protein[lookups.MASTER_INDEX]] = [protein]
            protein = next(all_protein_group_content)
        sorted_pgs = sort_protein_groups(proteins_in_groups)
	
        psm_masters = []
        psm_pg_proteins = []
        for master, group in sorted_pgs.items():
            psm_masters.append(master)
            psm_pg_proteins.append([protein[lookups.PROTEIN_ACC_INDEX] for protein in group])
        outpsm = {mzidtsvdata.HEADER_MASTER_PROT: ';'.join(psm_masters),
                  mzidtsvdata.HEADER_PG_CONTENT: ';'.join(
                  [','.join([y for y in x]) for x in psm_pg_proteins]),
                  mzidtsvdata.HEADER_PG_AMOUNT_PROTEIN_HITS: ';'.join(count_protein_group_hits(lineproteins, psm_pg_proteins))
                  }
        outpsm.update(psm)
        rownr += 1
        yield outpsm


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
    complete with sequences, psm_ids and scores. Master proteins are included
    in this group.
    
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
    If shared master proteins are found, report only the first,
    we will sort the whole proteingroup later anyway. In that
    case, the master reported here may be temporary."""
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


def count_protein_group_hits(lineproteins, groups):
    """Takes a list of protein accessions and a list of protein groups
    content from DB. Counts for each group in list how many proteins
    are found in lineproteins. Returns list of str amounts.
    """
    hits = []
    print('all groups:', groups)
    for group in groups:
        hits.append(0)
        print('line', lineproteins)
        for protein in lineproteins:
            if protein in group:
                hits[-1] += 1
    return [str(x) for x in hits]


def sort_protein_groups(pgroups):
    sortfnxs = [sort_pgroup_peptides,
                sort_pgroup_psms,
                sort_pgroup_score,
                sort_pgroup_coverage,
                sort_alphabet,
                ]
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
    return sort_amounts(proteins, lookups.PROTEIN_SCORE_INDEX)


def sort_pgroup_coverage(proteins):
    return [proteins]

def sort_evidence_score(proteins):
    return [proteins]

def sort_coin_toss(proteins):
    return [proteins]

def sort_alphabet(proteins):
    return [sorted(proteins, key=lambda x: x[2])]
