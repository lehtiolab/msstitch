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


def build_content_db(pgdb):
    use_evi = pgdb.check_evidence_tables()
    all_master_psm_proteins = pgdb.get_master_contentproteins_psms()
    all_master_psms = pgdb.get_all_master_psms()
    lastpsmmaster, masterpsm = next(all_master_psms)
    master_psms = {masterpsm}
    (lastcontentmaster, contentpsm, protein,
     pepseq, score, evid, cover) = next(all_master_psm_proteins)
    content = add_protein_psm_to_pre_proteingroup(dict(), protein, pepseq,
                                                  contentpsm, score, evid,
                                                  cover)
    protein_groups = []
    new_masters = {}
    for master, masterpsm in all_master_psms:
        # outer loop gets all master PSMs
        if master != lastpsmmaster:
            lastcontentmaster, pgroup, content = fetch_pg_content(
                all_master_psm_proteins, lastcontentmaster, lastpsmmaster,
                content, master_psms)
            new_master = sorters.sort_to_get_master(pgroup, use_evi)
            new_masters[new_master['master_id']] = new_master['protein_acc']
            pgroup = [[pg[2], pg[1], pg[3], pg[4], pg[5]] for pg in pgroup]
            protein_groups.extend(pgroup)
            master_psms = set()
            lastpsmmaster = master
        master_psms.add(masterpsm)
    lastcontentmaster, pgroup, content = fetch_pg_content(
        all_master_psm_proteins, lastcontentmaster,
        lastpsmmaster, content, master_psms)
    new_master = sorters.sort_to_get_master(pgroup, use_evi)
    new_masters[new_master['master_id']] = new_master['protein_acc']
    new_masters = ((acc, mid) for mid, acc in new_masters.items())
    pgroup = [[pg[2], pg[1], pg[3], pg[4], pg[5]] for pg in pgroup]
    protein_groups.extend(pgroup)
    pgdb.update_master_proteins(new_masters)
    pgdb.store_protein_group_content(protein_groups)
    pgdb.index_protein_group_content()


def fetch_pg_content(all_master_psm_proteins, lastcontentmaster, psmmaster,
                     content, master_psms):
    filtered = False
    for (contentmaster, contentpsm, protein, pepseq,
         score, evid, cover) in all_master_psm_proteins:
        # Inner loop gets protein group content from DB join table
        if contentmaster != lastcontentmaster:
            p_group = filter_proteins_with_missing_psms(content,
                                                        master_psms)
            lastcontentmaster, content = contentmaster, dict()
            content = add_protein_psm_to_pre_proteingroup(content,
                                                          protein,
                                                          pepseq,
                                                          contentpsm,
                                                          score,
                                                          evid,
                                                          cover)
            filtered = True
            break
        content = add_protein_psm_to_pre_proteingroup(content, protein,
                                                      pepseq,
                                                      contentpsm,
                                                      score, evid, cover)
    # The last protein-psm will not be caught by the if statement and thus
    # will there be one missing. But the break-construction makes that breaking
    # from the loop and loop exiting with StopIteration will both come to the
    # same point. Therefore we check if filtering has occurred before returning
    if not filtered:
        p_group = filter_proteins_with_missing_psms(content,
                                                    master_psms)
    pgroup_out = get_protein_group_content(p_group, psmmaster)
    return lastcontentmaster, pgroup_out, content


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


def filter_proteins_with_missing_psms(proteins, pg_psms):
    filtered_protein_map = {}
    for protein, protein_psms in proteins.items():
        filter_out = False
        for psm_id in [psm[0] for peptide in protein_psms.values()
                       for psm in peptide]:
            if psm_id not in pg_psms:
                filter_out = True
                break
        if filter_out:
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
