MZIDTSV_PEP_COL = 9
MZIDTSV_PROT_COL = 10
DB_STORE_CHUNK = 100000

from collections import OrderedDict

from app.readers import tsv as tsvreader
from app.actions.mzidtsv import proteingroup_sorters as sorters


def build_proteingroup_db(pgdb):
    build_master_db(pgdb)
    build_coverage(pgdb)
    build_content_db(pgdb)


def build_master_db(pgdb):
    psm_masters = OrderedDict()
    allmasters = {}
    allpsms = {}
    for psmid, protein in pgdb.get_all_psm_protein_relations():
        try:
            allpsms[psmid].append(protein)
        except KeyError:
            allpsms[psmid] = [protein]
    while len(allpsms) > 0:
        psm_id, proteins = allpsms.popitem()
        pepprotmap = pgdb.get_protpepmap_from_proteins(proteins)
        masters = get_masters(pepprotmap)
        for psm, master in [(p, m) for p, pmasters in masters.items()
                            for m in pmasters]:
            try:
                psm_masters[psm].add(master)
            except KeyError:
                psm_masters[psm] = set([master])
            allmasters[master] = 1
            if psm in allpsms:
                del(allpsms[psm])
    print('Collected {0} masters, {1} PSM-master mappings'.format(
        len(allmasters), len(psm_masters)))
    pgdb.store_masters(allmasters, psm_masters)


def process_pgroup_candidates(candidates, protein_psm_map):
    prepgroup = {}
    for candidate in candidates:
        master, psm_id, prot_id, seq, score, evid, cov = candidate
        prepgroup = add_protein_psm_to_pre_proteingroup(prepgroup, prot_id,
                                                        seq, psm_id, score,
                                                        evid, cov)
    pgroup = filter_proteins_with_missing_psms(prepgroup, protein_psm_map)
    return get_protein_group_content(pgroup, master)


def build_content_db(pgdb):
    protein_psms = {}
    for prot, psm in pgdb.get_protein_psm_records():
        try:
            protein_psms[prot].add(psm)
        except KeyError:
            protein_psms[prot] = set([psm])
    use_evi = pgdb.check_evidence_tables()
    pg_candidates = pgdb.get_protein_group_candidates()
    pre_protein_group = [next(pg_candidates)]
    lastmaster = pre_protein_group[0][0]
    protein_groups, new_masters = [], {}
    for protein_candidate in pg_candidates:
        if protein_candidate[0] != lastmaster:
            pgroup = process_pgroup_candidates(pre_protein_group, protein_psms)
            new_master = sorters.sort_to_get_master(pgroup, use_evi)
            new_masters[new_master['master_id']] = new_master['protein_acc']
            protein_groups.extend(pgroup)
            lastmaster, pre_protein_group = protein_candidate[0], []
        pre_protein_group.append(protein_candidate)
    pgroup = process_pgroup_candidates(pre_protein_group, protein_psms)
    new_master = sorters.sort_to_get_master(pgroup, use_evi)
    new_masters[new_master['master_id']] = new_master['protein_acc']
    protein_groups.extend(pgroup)
    protein_groups = [[pg[2], pg[1], pg[3], pg[4], pg[5]]
                      for pg in protein_groups]
    new_masters = ((acc, mid) for mid, acc in new_masters.items())
    pgdb.update_master_proteins(new_masters)
    pgdb.store_protein_group_content(protein_groups)
    pgdb.index_protein_group_content()


def add_protein_psm_to_pre_proteingroup(prepgmap, protein, pepseq,
                                        psm_id, score, evid, cover):
    protpsm_unit = (psm_id, float(score), evid, cover)
    try:
        prepgmap[protein][pepseq].add(protpsm_unit)
    except KeyError:
        try:
            prepgmap[protein][pepseq] = {protpsm_unit}
        except KeyError:
            prepgmap[protein] = {pepseq: {protpsm_unit}}
    return prepgmap


def filter_proteins_with_missing_psms(proteins, allprotein_psms):
    filtered_protein_map = {}
    for protein, protein_psms in proteins.items():
        protein_masterpsms = {psm[0] for peptide in protein_psms.values()
                              for psm in peptide}
        if allprotein_psms[protein].difference(protein_masterpsms):
            continue
        else:
            filtered_protein_map[protein] = protein_psms
    return filtered_protein_map


def build_coverage(pgdb):
    coverage = {}
    for acc, seq, psm_id, psmseq in pgdb.get_all_proteins_psms_seq():
        try:
            coverage[acc]['seq'] = seq
        except KeyError:
            coverage[acc] = {'seq': seq, 'psms': []}
        coverage[acc]['psms'].append(psmseq)
    pgdb.store_coverage(generate_coverage(coverage))


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
            if peps.issubset(subpeps):
                if peps.union(subpeps) > peps:
                    ismaster = False
                    break
                elif peps.intersection(subpeps) == peps:
                    multimaster.update({protein, subprotein})
        if not ismaster:
            continue
        elif multimaster:
            premaster = sorted(list(multimaster))[0]
        else:
            premaster = protein
        for pep in peps:
            try:
                masters[pep].add(premaster)
            except KeyError:
                masters[pep] = {premaster}
    return masters


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
            psmseq = tsvreader.strip_modifications(psmseq)
            # FIXME try block is for problems with coverage, see if it is
            # needed
            try:
                start = seq.index(psmseq)
            except:
                print('CANNOT FIND PSM seq {0} in seq {1} '
                      'for acc {2}'.format(psmseq, seq, acc))
            coverage_aa_indices.update(range(start, start + len(psmseq)))
        yield (acc, len(coverage_aa_indices) / len(seq))


def get_protein_group_content(pgmap, master):
    """For each master protein, we generate the protein group proteins
    complete with sequences, psm_ids and scores. Master proteins are included
    in this group.

    Returns a list of [protein, master, pep_hits, psm_hits, protein_score],
    which is ready to enter the DB table.
    """
    # first item (0) is only a placeholder so the lookup.INDEX things get the
    # correct number. Would be nice with a solution, but the INDEXes were
    # originally made for mzidtsv protein group adding.
    pg_content = [[0, master, protein, len(peptides), len([psm for pgpsms in
                                                           peptides.values()
                                                           for psm in pgpsms]),
                   sum([psm[1] for pgpsms in peptides.values()
                        for psm in pgpsms]),  # score
                   next(iter(next(iter(peptides.values()))))[3],  # coverage
                   next(iter(next(iter(peptides.values()))))[2],  # evid level
                   ]
                  for protein, peptides in pgmap.items()]
    return pg_content
