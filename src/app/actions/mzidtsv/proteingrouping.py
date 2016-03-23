from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata
from app.lookups.sqlite import proteingroups as lookups
from app.actions.mzidtsv import proteingroup_sorters as sorters


def get_header_with_proteingroups(header):
    ix = header.index(mzidtsvdata.HEADER_PROTEIN) + 1
    return header[:ix] + mzidtsvdata.HEADER_PG + header[ix:]


def get_all_proteins_from_unrolled_psm(rownr, pgdb):
    return pgdb.get_proteins_for_peptide([rownr])


def generate_psms_with_proteingroups(fn, oldheader, newheader, pgdb,
                                     unroll=False):
    rownr = 0
    use_evi = pgdb.check_evidence_tables()
    all_protein_group_content = pgdb.get_all_psms_proteingroups(use_evi)
    protein = next(all_protein_group_content)
    for psm in tsvreader.generate_tsv_psms(fn, oldheader):
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
        sorted_pgs = sorters.sort_protein_groups(proteins_in_groups, use_evi)
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
