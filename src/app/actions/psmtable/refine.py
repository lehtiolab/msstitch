import re
from collections import OrderedDict

from app.readers import tsv as tsvreader
from app.readers import fasta as fastareader
from app.dataformats import mzidtsv as mzidtsvdata

from app.lookups.sqlite import psms as lookups

DB_STORE_CHUNK = 100000


def create_header(oldheader, genes, proteingroup, precursor, isob_header, bioset,
        miscleav, specfncolnr):
    header = oldheader[:]
    p_ix = header.index(mzidtsvdata.HEADER_PROTEIN) + 1
    if genes:
        newfields = [mzidtsvdata.HEADER_GENE, mzidtsvdata.HEADER_SYMBOL,
                     mzidtsvdata.HEADER_DESCRIPTION]
        header = header[:p_ix] + newfields + oldheader[p_ix:]
    if proteingroup:
        header = header[:p_ix] + mzidtsvdata.HEADER_PG + header[p_ix:]
    if precursor:
        header += [mzidtsvdata.HEADER_PRECURSOR_QUANT]
    if isob_header:
        header += isob_header
    psmdatafields = mzidtsvdata.MOREDATA_HEADER
    if bioset:
        psmdatafields = [mzidtsvdata.HEADER_SETNAME] + psmdatafields
    if miscleav:
        psmdatafields.append(mzidtsvdata.HEADER_MISSED_CLEAVAGE)
    header = header[:specfncolnr +1] + psmdatafields + header[specfncolnr + 1:]
    return header


def create_psm_lookup(fn, fastafn, header, pgdb, unroll, specfncol, fastadelim, genefield):
    """Reads PSMs from file, stores them to a database backend in chunked PSMs.
    """
    proteins = store_proteins_descriptions(pgdb, fastafn, fn, header,
                                           fastadelim, genefield)
    mzmlmap = pgdb.get_mzmlfile_map()
    sequences = {}
    for psm in tsvreader.generate_tsv_psms(fn, header):
        seq = tsvreader.get_psm_sequence(psm, unroll)
        sequences[seq] = 1
    pgdb.store_pepseqs(((seq,) for seq in sequences))
    pepseqmap = pgdb.get_peptide_seq_map()
    psms = []
    for row, psm in enumerate(tsvreader.generate_tsv_psms(fn, header)):
        specfn, psm_id, specscanid, seq, score = tsvreader.get_psm(psm, unroll, specfncol)
        if len(psms) % DB_STORE_CHUNK == 0:
            pgdb.store_psms(psms)
            psms = []
        psms.append({'rownr': row,
                     'psm_id': psm_id,
                     'seq': pepseqmap[seq],
                     'score': score,
                     'specfn': mzmlmap[specfn],
                     'spec_id': '{}_{}'.format(mzmlmap[specfn], specscanid),
                     })
    pgdb.store_psms(psms)
    pgdb.index_psms()
    store_psm_protein_relations(fn, header, pgdb, proteins)


def store_proteins_descriptions(pgdb, fastafn, tsvfn, header, fastadelim, genefield):
    if not fastafn:
        prots = {}
        for psm in tsvreader.generate_tsv_psms(tsvfn, header):
            prots.update({x: 1 for x in
                             tsvreader.get_proteins_from_psm(psm)})
        prots = [(protein,) for protein in prots.keys()]
        pgdb.store_proteins(prots)
    else:
        prots, seqs, desc, evids, ensgs, symbols = fastareader.get_proteins_for_db(
            fastafn, fastadelim, genefield)
        pgdb.store_fasta(prots, evids, seqs, desc, ensgs, symbols)
    return set([x[0] for x in prots])


def store_psm_protein_relations(fn, header, pgdb, proteins):
    """Reads PSMs from file, extracts their proteins and peptides and passes
    them to a database backend in chunks.
    """
    # TODO do we need an OrderedDict or is regular dict enough?
    # Sorting for psm_id useful?
    allpsms = OrderedDict()
    last_id, psmids_to_store = None, set()
    store_soon = False
    for psm in tsvreader.generate_tsv_psms(fn, header):
        psm_id, prots = tsvreader.get_pepproteins(psm)
        # TODO can this be removed permanently? 
        # Filter proteins to only include those that match the protein 
        # accessions in fasta so we get the correct names, filter out the badly annotated peptides
        # prots = [x for x in prots if x in proteins]
        try:
            # In case the PSMs are presented unrolled
            allpsms[psm_id].extend(prots)
        except KeyError:
            allpsms[psm_id] = prots
        if len(psmids_to_store) > DB_STORE_CHUNK:
            store_soon = True
        if store_soon and last_id != psm_id:
            pgdb.store_peptides_proteins(allpsms, psmids_to_store)
            store_soon = False
            psmids_to_store = set()
        psmids_to_store.add(psm_id)
        last_id = psm_id
    if len(psmids_to_store) > 0:
        pgdb.store_peptides_proteins(allpsms, psmids_to_store)
    pgdb.index_protein_peptides()
    return allpsms

def add_genes_to_psm_table(psms, pgdb):
    gpmap = pgdb.get_protein_gene_map()
    for psm in psms:
        outpsm = {x: y for x, y in psm.items()}
        proteins = tsvreader.get_proteins_from_psm(psm)
        outpsm[mzidtsvdata.HEADER_GENE] = ';'.join(get_genes(proteins, gpmap))
        symbols = get_symbols(proteins, gpmap)
        desc = get_descriptions(proteins, gpmap)
        outpsm[mzidtsvdata.HEADER_SYMBOL] = ';'.join(symbols)
        outpsm[mzidtsvdata.HEADER_DESCRIPTION] = ';'.join(desc)
        yield outpsm


def get_mapped(proteins, gpmap, outtype):
    # FIXME multiple vals for a protein?
    outvals = [gpmap[protein][outtype] for protein in proteins]
    outvals = OrderedDict([(val, 1) for val in outvals]).keys()
    if None in outvals:
        return ['NA']
    else:
        return outvals


def get_genes(proteins, gpmap):
    return get_mapped(proteins, gpmap, 'gene')


def get_descriptions(proteins, gpmap):
    descriptions = get_mapped(proteins, gpmap, 'desc')
    return [x.replace('\t', ' ') for x in descriptions]


def get_symbols(proteins, gpmap):
    return get_mapped(proteins, gpmap, 'symbol')


def get_all_proteins_from_unrolled_psm(rownr, pgdb):
    return pgdb.get_proteins_for_peptide([rownr])


# FIXME tsv
def generate_psms_with_proteingroups(psms, pgdb, unroll=False):
    rownr = 0
    use_evi = pgdb.check_evidence_tables()
    all_protein_group_content = pgdb.get_all_psms_proteingroups(use_evi)
    protein = next(all_protein_group_content)
    for psm in psms:
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
        sorted_pgs = sort_protein_groups(proteins_in_groups, use_evi)
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


def get_sortfnxs(evidence):
    sortfnxs = [sort_pgroup_peptides,
                sort_pgroup_psms,
                sort_pgroup_score,
                ]
    if evidence:
        sortfnxs.append(sort_evidence_score)
    sortfnxs.extend([sort_pgroup_coverage, sort_alphabet])
    return sortfnxs


def sort_to_get_master(pgroup, evidence):
    sortfnxs = get_sortfnxs(evidence)
    sorted_pg = sort_protein_group(pgroup, sortfnxs, 0)
    return {'master_id': sorted_pg[0][lookups.MASTER_INDEX],
            'protein_acc': sorted_pg[0][lookups.PROTEIN_ACC_INDEX]}


def sort_protein_groups(pgroups, evidence):
    """Gets a protein groups containing dict pgroups, for each master (key)
    there is a list of protein tuples that comprise a proteingroup. This loops
    the groups and returns a sorted group. Assumes master is already sorted
    as the proper master"""
    sortfnxs = get_sortfnxs(evidence)
    pgroups_out = {}
    for master, pgroup in pgroups.items():
        sorted_pg = sort_protein_group(pgroup, sortfnxs, 0)
        pgroups_out[master] = sorted_pg
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


def sort_amounts(proteins, sort_index):
    """Generic function for sorting peptides and psms. Assumes a higher
    number is better for what is passed at sort_index position in protein."""
    amounts = {}
    for protein in proteins:
        amount_x_for_protein = protein[sort_index]
        try:
            amounts[amount_x_for_protein].append(protein)
        except KeyError:
            amounts[amount_x_for_protein] = [protein]
    return [v for k, v in sorted(amounts.items(), reverse=True)]


def sort_pgroup_peptides(proteins):
    return sort_amounts(proteins, lookups.PEPTIDE_COUNT_INDEX)


def sort_pgroup_psms(proteins):
    return sort_amounts(proteins, lookups.PSM_COUNT_INDEX)


def sort_pgroup_score(proteins):
    return sort_amounts(proteins, lookups.PROTEIN_SCORE_INDEX)


def sort_pgroup_coverage(proteins):
    return sort_amounts(proteins, lookups.COVERAGE_INDEX)


def sort_evidence_score(proteins):
    return sort_amounts(proteins, lookups.EVIDENCE_LVL_INDEX)


def sort_alphabet(proteins):
    return [sorted(proteins, key=lambda x: x[2])]


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
            new_master = sort_to_get_master(pgroup, use_evi)
            new_masters[new_master['master_id']] = new_master['protein_acc']
            protein_groups.extend(pgroup)
            lastmaster, pre_protein_group = protein_candidate[0], []
        pre_protein_group.append(protein_candidate)
    pgroup = process_pgroup_candidates(pre_protein_group, protein_psms)
    new_master = sort_to_get_master(pgroup, use_evi)
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


def generate_psms_quanted(quantdb, psms, isob_header, isobaric=False, precursor=False):
    """Takes dbfn and connects, gets quants for each line in tsvfn, sorts
    them in line by using keys in quantheader list."""
    allquants, sqlfields = quantdb.select_all_psm_quants(isobaric, precursor)
    quant = next(allquants)
    for rownr, psm in enumerate(psms):
        outpsm = {x: y for x, y in psm.items()}
        if precursor:
            pquant = quant[sqlfields['precursor']]
            if pquant is None:
                pquant = 'NA'
            outpsm.update({mzidtsvdata.HEADER_PRECURSOR_QUANT: str(pquant)})
        if isobaric:
            isoquants = {}
            while quant[0] == rownr:
                isoquants.update({quant[sqlfields['isochan']]:
                                  str(quant[sqlfields['isoquant']])})
                try:
                    quant = next(allquants)
                except StopIteration:
                    # last PSM, break from while loop or it is not yielded at all
                    break
            outpsm.update(get_quant_NAs(isoquants, isob_header))
        else:
            try:
                quant = next(allquants)
            except StopIteration:
                # last PSM, needs explicit yield/break or it will not be yielded
                yield outpsm
                break
        yield outpsm


def get_quant_NAs(quantdata, quantheader):
    """Takes quantdata in a dict and header with quantkeys
    (eg iTRAQ isotopes). Returns dict of quant intensities
    with missing keys set to NA."""
    out = {}
    for qkey in quantheader:
        out[qkey] = quantdata.get(qkey, 'NA')
    return out


def count_missed_cleavage(full_pepseq, count=0):
    '''Regex .*[KR][^P] matches until the end and checks if there is a final
    charachter so this will not match the tryptic residue'''
    pepseq = re.sub('[\+\-]\d*.\d*', '', full_pepseq)
    match = re.match('.*[KR][^P]', pepseq)
    if match:
        count += 1
        return count_missed_cleavage(match.group()[:-1], count)
    else:
        return count


def generate_psms_spectradata(lookup, psms, bioset, miscleav):
    psm_specdata = zip(enumerate(psms), lookup.get_exp_spectra_data_rows())
    for (row, psm), specdata in psm_specdata:
        outpsm = {x: y for x, y in psm.items()}
        if row == int(specdata[0]):
            inj = str(specdata[3])
            imob = str(specdata[4])
            outpsm.update({mzidtsvdata.HEADER_RETENTION_TIME: str(specdata[2]),
                           mzidtsvdata.HEADER_INJECTION_TIME: inj if inj != 'None' else 'NA',
                           mzidtsvdata.HEADER_ION_MOB: imob if imob != 'None' else 'NA',
                           })
            if bioset:
                outpsm[mzidtsvdata.HEADER_SETNAME] = str(specdata[1])
        else:
            raise RuntimeError('PSM with row nr {} has no rownr in DB. '
                               'Current DB row is {}'.format(row, specdata[0]))
        if miscleav:
            outpsm[mzidtsvdata.HEADER_MISSED_CLEAVAGE] = count_missed_cleavage(outpsm[mzidtsvdata.HEADER_PEPTIDE])
        yield outpsm
