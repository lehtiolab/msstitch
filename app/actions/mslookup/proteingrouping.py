MZIDTSV_PEP_COL = 9
MZIDTSV_PROT_COL = 10
DB_STORE_CHUNK = 500000

from collections import OrderedDict

from app.lookups.sqlite.proteingroups import ProteinGroupDB
from app.readers import tsv as tsvreader
from app.readers import fasta as fastareader
from app.actions.mzidtsv import confidencefilters as conffilt


def create_protein_pep_lookup(fn, header, pgdb, confkey, conflvl,
                              lower_is_better, unroll=False, fastafn=None,
                              specfncol=None):
    """Reads PSMs from file, extracts their proteins and peptides and passes
    them to a database backend in chunked PSMs.
    """
    mzmlmap = pgdb.get_mzmlfile_map()
    proteins, sequences, evidences = fastareader.get_proteins_for_db(
        fastafn)
    pgdb.store_proteins(proteins, evidences, sequences)
    protein_descriptions = fastareader.get_proteins_descriptions(fastafn)
    pgdb.store_descriptions(protein_descriptions)
    
    rownr, last_id, peptides_proteins = 0, None, {}
    store_soon = False
    for psm in tsvreader.generate_tsv_lines_multifile(fn, header):
        if not conffilt.passes_filter(psm, conflvl, confkey, lower_is_better):
            rownr += 1
            continue
        specfn, scan, seq, score, prots = tsvreader.get_pepproteins(
            psm, unroll, specfncol)
        psm_id = tsvreader.get_psm_id(psm)
        if peptides_proteins and len(peptides_proteins) % DB_STORE_CHUNK == 0:
            store_soon = True
        if store_soon and last_id != psm_id:
            pgdb.store_peptides_proteins(peptides_proteins)
            store_soon = False
            peptides_proteins = {}
        peptides_proteins[rownr] = {'psm_id': psm_id,
                                    'seq': seq,
                                    'proteins': prots,
                                    'score': score,
                                    'specfn': mzmlmap[specfn],
                                    'scannr': scan,
                                    }
        last_id = psm_id
        rownr += 1
    pgdb.store_peptides_proteins(peptides_proteins)
    pgdb.index_protein_peptides()
    return pgdb


def build_proteingroup_db(fn, oldheader, pgdb, confkey, conflvl,
                          lower_is_better, unroll, coverage):
    build_master_db(fn, oldheader, pgdb, confkey, conflvl, lower_is_better,
                    unroll)
    build_content_db(pgdb)
    if coverage:
        build_coverage(pgdb)


def build_master_db(fn, oldheader, pgdb, confkey, conflvl, lower_is_better,
                    unroll):
    psm_masters = OrderedDict()
    allmasters = {}
    rownr = 0
    for line in tsvreader.generate_tsv_lines_multifile(fn, oldheader):
        if not conffilt.passes_filter(line, conflvl, confkey, lower_is_better):
            rownr += 1
            continue
        psm_id = tsvreader.get_psm_id(line)
        if unroll:
            lineproteins = pgdb.get_proteins_for_peptide([psm_id])
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
