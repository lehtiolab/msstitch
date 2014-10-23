from collections import OrderedDict

from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata
import app.sqlite as lookups
import app.preparation.mzidtsv.proteingroup_sorters as sorters
from app.preparation.mzidtsv import confidencefilters as conffilt


def get_header_with_proteingroups(header):
    ix = header.index(mzidtsvdata.HEADER_PROTEIN) + 1
    return header[:ix] + mzidtsvdata.HEADER_PG + header[ix:]


def get_all_proteins_from_unrolled_psm(rownr, pgdb):
    return pgdb.get_proteins_for_peptide([rownr])


def generate_psms_with_proteingroups(fn, oldheader, newheader, pgdb, confkey,
                                     conflvl, lower_is_better, unroll=False,
                                     coverage=False, evidence_levels=False):
    build_master_db(fn, oldheader, pgdb, confkey, conflvl, lower_is_better,
                    unroll)
    build_content_db(pgdb)
    if coverage:
        build_coverage(pgdb)
    rownr = 0
    all_protein_group_content = pgdb.get_all_psms_proteingroups(
        coverage, evidence_levels)
    protein = next(all_protein_group_content)
    for psm in tsvreader.generate_tsv_psms(fn, oldheader):
        #if rownr % 10000 == 0 and rownr != 0:
        #    break
        if not conffilt.passes_filter(psm, conflvl, confkey, lower_is_better):
            rownr += 1
            continue
        if unroll:
            psm_id = tsvreader.get_psm_id(psm)
            lineproteins = get_all_proteins_from_unrolled_psm(psm_id, pgdb)
        else:
            lineproteins = tsvreader.get_proteins_from_psm(psm)
        proteins_in_groups = {}
        while protein[0] == rownr:
            try:
                proteins_in_groups[protein[
                    lookups.MASTER_INDEX]].append(protein)
            except KeyError:
                proteins_in_groups[protein[lookups.MASTER_INDEX]] = [protein]
            try:
                protein = next(all_protein_group_content)
            except StopIteration:
                protein = [-1]
        sorted_pgs = sorters.sort_protein_groups(proteins_in_groups, coverage,
                                                 evidence_levels)
        psm_masters = []
        psm_pg_proteins = []
        for master, group in sorted_pgs.items():
            psm_masters.append(master)
            psm_pg_proteins.append([protein[lookups.PROTEIN_ACC_INDEX]
                                    for protein in group])
        outpsm = {mzidtsvdata.HEADER_MASTER_PROT: ';'.join(psm_masters),
                  mzidtsvdata.HEADER_PG_CONTENT: ';'.join(
                      [','.join([y for y in x]) for x in psm_pg_proteins]),
                  mzidtsvdata.HEADER_PG_AMOUNT_PROTEIN_HITS: ';'.join(
                      count_protein_group_hits(lineproteins, psm_pg_proteins))
                  }
        outpsm.update(psm)
        rownr += 1
        yield outpsm


def build_master_db(fn, oldheader, pgdb, confkey, conflvl, lower_is_better,
                    unroll):
    psm_masters = OrderedDict()
    allmasters = {}
    rownr = 0
    for line in tsvreader.generate_tsv_psms(fn, oldheader):
        if not conffilt.passes_filter(line, conflvl, confkey, lower_is_better):
            rownr += 1
            continue
        psm_id = tsvreader.get_psm_id(line)
        if unroll:
            lineproteins = get_all_proteins_from_unrolled_psm(psm_id, pgdb)
        else:
            lineproteins = tsvreader.get_proteins_from_psm(line)
        pepprotmap = pgdb.get_protpepmap_from_proteins(lineproteins)
        masters = get_masters(pepprotmap)
        psm_masters[psm_id] = {x: 1 for x in masters}
        allmasters.update({x: 1 for x in masters})
        rownr += 1
    pgdb.store_masters(allmasters, psm_masters)
    return allmasters


def build_content_db(pgdb):
    protein_groups = []
    allpsms_masters = pgdb.get_allpsms_masters()
    lastmaster, psms = next(allpsms_masters)
    psms = [psms]
    for master, psm in allpsms_masters:
        if master == lastmaster:
            psms.append(psm)
        else:
            protein_groups.extend(get_protein_group_content(lastmaster, psms,
                                                            pgdb))
            psms = [psm]
        lastmaster = master
    pgdb.store_protein_group_content(protein_groups)


def build_coverage(pgdb):
    coverage = {}
    for acc, seq, psm_id, psmseq in pgdb.get_all_proteins_psms_seq():
        try:
            coverage[acc]['seq'] = seq
        except KeyError:
            coverage[acc] = {'seq': seq, 'psms': []}
        coverage[acc]['psms'].append(psmseq)
    pgdb.store_coverage(generate_coverage(coverage))


def generate_coverage(seqinfo):
    """From a dict containing protein accessions and sequences/PSM sequences,
    this function returns a generator that calculates coverages for each
    protein and returns the accession and coverage percentage.
    Coverage is done by finding peptides in the protein seq using seq.index
    and marking the range. May be slow."""
    for acc, protinfo in seqinfo.items():
        coverage_aa_indices = set()
        seq = protinfo['seq']
        for psmseq in protinfo['psms']:
            start = seq.index(psmseq)
            coverage_aa_indices.update(range(start, start + len(psmseq)))
        yield (acc, len(coverage_aa_indices) / len(seq))


def get_protein_group_content(master, psms, pgdb):
    """For each master protein, we generate the protein group proteins
    complete with sequences, psm_ids and scores. Master proteins are included
    in this group.

    Returns a list of [protein, master, pep_hits, psm_hits, protein_score],
    which is ready to enter the DB table.
    """
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
    pgmap = [[protein, master, len(peptides), len([psm for pgpsms in
                                                   peptides.values()
                                                   for psm in pgpsms]),
              sum([psm[1] for pgpsms in peptides.values() for psm in pgpsms])]
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
    for group in groups:
        hits.append(0)
        for protein in lineproteins:
            if protein in group:
                hits[-1] += 1
    return [str(x) for x in hits]
