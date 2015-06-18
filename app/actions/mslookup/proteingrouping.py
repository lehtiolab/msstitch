MZIDTSV_PEP_COL = 9
MZIDTSV_PROT_COL = 10
DB_STORE_CHUNK = 100000

from collections import OrderedDict

from app.lookups.sqlite.proteingroups import ProteinGroupDB
from app.readers import tsv as tsvreader
from app.readers import fasta as fastareader
from app.actions.mzidtsv import confidencefilters as conffilt


def create_protein_pep_lookup(fn, header, pgdb, confkey, conflvl,
                              lower_is_better, fastafn=None):
    """Reads PSMs from file, extracts their proteins and peptides and passes
    them to a database backend in chunks.
    """
    proteins, sequences, evidences = fastareader.get_proteins_for_db(
        fastafn)
    pgdb.store_proteins(proteins, evidences, sequences)
    protein_descriptions = fastareader.get_proteins_descriptions(fastafn)
    pgdb.store_descriptions(protein_descriptions)
    # TODO do we need an OrderedDict or is regular dict enough? Sorting for psm_id useful? 
    allpsms = OrderedDict()
    last_id, psmids_to_store = None, set()
    store_soon = False
    for psm in tsvreader.generate_tsv_lines_multifile(fn, header):
        if not conffilt.passes_filter(psm, conflvl, confkey, lower_is_better):
            continue
        psm_id, prots = tsvreader.get_pepproteins(psm)
        try:
            allpsms[psm_id].extend(prots)
        except KeyError:
            allpsms[psm_id] = prots
        if len(psmids_to_store) % DB_STORE_CHUNK == 0:
            store_soon = True
        if store_soon and last_id != psm_id:
            pgdb.store_peptides_proteins(allpsms, psmids_to_store)
            store_soon = False
            psmids_to_store = set()
        psmids_to_store.add(psm_id)
        last_id = psm_id
    pgdb.store_peptides_proteins(allpsms, psmids_to_store)
    pgdb.index_protein_peptides()
    return allpsms


def build_proteingroup_db(pgdb, allpsms,
                          coverage):
    build_master_db(pgdb, allpsms)
    build_content_db(pgdb)
    if coverage:
        build_coverage(pgdb)


def build_master_db(pgdb, allpsms):
    psm_masters = OrderedDict()
    allmasters = {}
    while len(allpsms) > 0:
        psm_id, proteins = allpsms.popitem()
        pepprotmap = pgdb.get_protpepmap_from_proteins(proteins)
        masters = get_masters(pepprotmap)
        for psm, master in [(p, m) for p, masters in masters.items() for m in masters]:
            try:
                psm_masters[psm].add(master)
            except KeyError:
                psm_masters[psm] = set([master])
            allmasters[master] = 1
            if psm in allpsms:
                del(allpsms[psm])
    print('Collected {0} masters, {1} PSM-master mappings'.format(len(allmasters), len(psm_masters)))
    pgdb.store_masters(allmasters, psm_masters)


def build_content_db(pgdb):
    all_master_psm_proteins = pgdb.get_master_contentproteins_psms()
    all_master_psms = pgdb.get_all_master_psms()
    lastpsmmaster, masterpsm = next(all_master_psms)
    lastcontentmaster, contentpsm, protein, pepseq, score = next(all_master_psm_proteins)
    master_psms = {masterpsm}
    contentmap = add_protein_psm_to_pre_proteingroup(dict(), protein, pepseq, contentpsm, score)
    protein_groups = []
    for master, masterpsm in all_master_psms:
        # outer loop gets all master PSMs
        if master != lastpsmmaster:
            for contentmaster, contentpsm, protein, pepseq, score in all_master_psm_proteins:
                # Inner loop gets protein group content from another DB join table
                if contentmaster != lastcontentmaster:
                    proteingroup = filter_proteins_with_missing_psms(contentmap, master_psms)
                    lastcontentmaster, contentmap = contentmaster, dict()
                    contentmap = add_protein_psm_to_pre_proteingroup(contentmap, protein, pepseq, contentpsm, score)
                    break 
                contentmap = add_protein_psm_to_pre_proteingroup(contentmap, protein, pepseq, contentpsm, score)
            protein_groups.extend(get_protein_group_content(proteingroup, lastpsmmaster))
            master_psms = set()
            lastpsmmaster = master
        master_psms.add(masterpsm)
    pgdb.store_protein_group_content(protein_groups)


def add_protein_psm_to_pre_proteingroup(prepgmap, protein, pepseq, psm_id, score):
    score = float(score)
    try:
        prepgmap[protein][pepseq].add((psm_id, score))
    except KeyError:
        try:
            prepgmap[protein][pepseq] = {(psm_id, score)}
        except KeyError:
            prepgmap[protein] = {pepseq: {(psm_id, score)}}
    return prepgmap
                

def filter_proteins_with_missing_psms(proteins, pg_psms):
    filtered_protein_map = {}
    for protein, protein_psms in proteins.items():
        filter_out = False
        for psm_id, score in [psm for peptide in protein_psms.values() for psm in peptide]:
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
                print('CANNOT FIND PSM seq {0} in seq {1} for acc {2}'.format(psmseq, seq, acc))
            coverage_aa_indices.update(range(start, start + len(psmseq)))
        yield (acc, len(coverage_aa_indices) / len(seq))


def get_protein_group_content(pgmap, master):
    """For each master protein, we generate the protein group proteins
    complete with sequences, psm_ids and scores. Master proteins are included
    in this group.

    Returns a list of [protein, master, pep_hits, psm_hits, protein_score],
    which is ready to enter the DB table.
    """
    pgmap = [[protein, master, len(peptides), len([psm for pgpsms in
                                                   peptides.values()
                                                   for psm in pgpsms]),
              sum([psm[1] for pgpsms in peptides.values() for psm in pgpsms])]
             for protein, peptides in pgmap.items()]
    return pgmap
