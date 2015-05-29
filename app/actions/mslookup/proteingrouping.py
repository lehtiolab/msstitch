MZIDTSV_PEP_COL = 9
MZIDTSV_PROT_COL = 10
DB_STORE_CHUNK = 500000

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
    build_new_master_db(pgdb, allpsms)
    build_content_db(pgdb)
    if coverage:
        build_coverage(pgdb)

def build_new_master_db(pgdb, allpsms):
    psm_masters = OrderedDict()
    allmasters = {}
    while len(allpsms) > 0:
        psm_id, proteins = allpsms.popitem()
        pepprotmap = pgdb.get_protpepmap_from_proteins(proteins)
        masters = get_masters(pepprotmap)
        for psm_id_delete in [x for y in pepprotmap.values() for x in y]:
            try:
                psm_masters[psm_id_delete].update({x: 1 for x in masters})
            except KeyError:
                psm_masters[psm_id_delete] = {x: 1 for x in masters}
            if psm_id_delete in allpsms:
                del(allpsms[psm_id_delete])
        psm_masters[psm_id] = {x: 1 for x in masters}
        allmasters.update({x: 1 for x in masters})
    print('Collected {0} masters, {1} PSM-master mappings'.format(len(allmasters), len(psm_masters)))
    pgdb.store_masters(allmasters, psm_masters)

def build_master_db(pgdb):
    psm_masters = OrderedDict()
    allmasters = {}
    allpepprots = pgdb.get_all_pepprots()
    current_psm, protein_acc, prot_psm_id = next(allpepprots)
    pepprotmap = {protein_acc: [prot_psm_id]}
    for psm_id, protein_acc, prot_psm_id in allpepprots:
        if psm_id != current_psm:
            # got all psm_protein mappings for this psm
            # flush mappings and store masters, psm_master mapping in variables
            masters = get_masters(pepprotmap)
            psm_masters[current_psm] = {x: 1 for x in masters}
            allmasters.update({x: 1 for x in masters})
            pepprotmap = {}
        current_psm = psm_id
        try:
            pepprotmap[protein_acc].append(prot_psm_id)
        except KeyError:
            pepprotmap[protein_acc] = [prot_psm_id]
    print('Collected {0} masters, {1} PSM-master mappings'.format(len(allmasters), len(psm_masters)))
    pgdb.store_masters(allmasters, psm_masters)


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
    protein_groups.extend(get_protein_group_content(lastmaster, psms, pgdb))
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
