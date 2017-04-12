from collections import OrderedDict


from app.dataformats import peptable as peptabledata
from app.actions.mergetable import (simple_val_fetch, fill_mergefeature,
                                    get_isobaric_quant)
from app.actions.proteindata import create_featuredata_map


def build_peptidetable(pqdb, headerfields, isobaric=False,
                       precursor=False, fdr=False, pep=False,
                       genecentric=False):
    """Fetches peptides and quants from joined lookup table, loops through
    them and when all of a peptides quants/data have been collected, yields
    peptide quant information."""
    peptidedatamap = create_featuredata_map(pqdb, genecentric=genecentric,
                                            fill_fun=add_record_to_peptidedata)
    count_psms(peptidedatamap)
    empty_return = lambda x, y, z: {}
    iso_fun = {True: get_isobaric_quant, False: empty_return}[isobaric]
    ms1_fun = {True: get_precursor_quant, False: empty_return}[precursor]
    fdr_fun = {True: get_pep_fdr,
               False: empty_return}[fdr]
    pep_fun = {True: get_peptide_pep,
               False: empty_return}[pep]
    pdata_fun = get_protein_data
    peptide_sql, sqlfieldmap = pqdb.prepare_mergetable_sql(precursor, isobaric,
                                                           probability=False,
                                                           fdr=fdr, pep=pep)
    peptides = pqdb.get_merged_features(peptide_sql)
    peptide = next(peptides)
    outpeptide = {peptabledata.HEADER_PEPTIDE: peptide[sqlfieldmap['p_acc']]}
    check_pep = {k: v for k, v in outpeptide.items()}
    fill_mergefeature(outpeptide, iso_fun, ms1_fun, empty_return, fdr_fun,
                      pep_fun, pdata_fun, peptide, sqlfieldmap,
                      headerfields, peptidedatamap)
    for peptide in peptides:
        p_seq = peptide[sqlfieldmap['p_acc']]
        if p_seq != outpeptide[peptabledata.HEADER_PEPTIDE]:
            if outpeptide != check_pep:
                yield outpeptide
            outpeptide = {peptabledata.HEADER_PEPTIDE: p_seq}
            check_pep = {k: v for k, v in outpeptide.items()}
        fill_mergefeature(outpeptide, iso_fun, ms1_fun, empty_return, fdr_fun,
                          pep_fun, pdata_fun, peptide, sqlfieldmap,
                          headerfields, peptidedatamap)
    yield outpeptide


def get_precursor_quant(peptide, sqlmap, headerfields):
    return simple_val_fetch(peptide, sqlmap,
                            headerfields['precursorquant'][
                                peptabledata.HEADER_AREA], 'preq_val')


def get_pep_fdr(peptide, sqlmap, headerfields):
    return simple_val_fetch(peptide, sqlmap,
                            headerfields['peptidefdr'][
                                peptabledata.HEADER_QVAL], 'fdr_val')


def get_peptide_pep(peptide, sqlmap, headerfields):
    return simple_val_fetch(peptide, sqlmap,
                            headerfields['peptidepep'][
                                peptabledata.HEADER_PEP], 'pep_val')


def get_no_psms(peptide, pdata, headerfields):
    hfields = [peptabledata.HEADER_NO_PSM,
               ]
    seq = peptide[peptabledata.HEADER_PEPTIDE]
    outdict = {}
    for pool, psms in pdata[seq]['psms'].items():
        pool_values = [psms]
        outdict.update({headerfields['nopsms'][hfield][pool]: val
                        for (hfield, val) in zip(hfields, pool_values)})
    return outdict


def get_protein_data(peptide, pdata, headerfields):
    """These fields are currently not pool dependent so headerfields
    is ignored"""
    report = get_proteins(peptide, pdata, headerfields)
    return get_cov_descriptions(peptide, pdata, report)


def get_proteins(peptide, pdata, headerfields):
    seq = peptide[peptabledata.HEADER_PEPTIDE]
    outdict = {}
    try:
        proteins = ';'.join([x[0] for x in pdata[seq]['proteins']])
    except TypeError:
        pass
    else:
        outdict = {peptabledata.HEADER_PROTEINS: proteins}
    for pool, psms in pdata[seq]['psms'].items():
        outdict.update(
            {headerfields['proteindata'][peptabledata.HEADER_NO_PSM][pool]:
             psms})
    return outdict


def get_cov_descriptions(peptide, pdata, report):
    def format_val(value, valtype):
        formatter = {int: lambda x: str(x),
                     float: lambda x: str(float(x)),
                     str: lambda x: x,
                     }
        return formatter[valtype](value)

    seq = peptide[peptabledata.HEADER_PEPTIDE]
    for idx, key, keytype in zip([1, 2, 3, 4, 5],
                                 [peptabledata.HEADER_COVERAGES,
                                  peptabledata.HEADER_DESCRIPTIONS,
                                  peptabledata.HEADER_GENES,
                                  peptabledata.HEADER_ASSOCIATED,
                                  peptabledata.HEADER_NO_CONTENTPROTEINS],
                                 [float, str, str, str, int]):
        try:
            replist = [format_val(x[idx], keytype)
                       for x in pdata[seq]['proteins']]
        except (TypeError, IndexError):
            # index too high, or None skip, item will be NA'ed downstream
            continue
        if key in [peptabledata.HEADER_GENES, peptabledata.HEADER_ASSOCIATED]:
            replist = OrderedDict([(x, 1) for x in replist]).keys()
        try:
            report[key] = ';'.join(replist)
        except TypeError:
            # None in list, skip, NA will be output
            continue
    return report


def count_psms(pdata):
    for seq in pdata:
        for pool in pdata[seq]['psms']:
            pdata[seq]['psms'][pool] = len(pdata[seq]['psms'][pool])


def add_record_to_peptidedata(peptidedata, p_acc, pool, psmdata, genecentric,
                              pgcontentmap=None):
    seq, psm_id = psmdata[2], psmdata[3]
    if genecentric == 'plain':
        gene, desc, assoc_id, cov, pgcontent = None, None, None, None, None
    elif genecentric:
        gene = p_acc
        desc, assoc_id = psmdata[4], psmdata[5]
        cov, pgcontent = None, None
    else:
        desc, cov = psmdata[4], psmdata[5]
        gene, assoc_id = psmdata[6], psmdata[7]
        pgcontent = pgcontentmap[p_acc]
    try:
        peptidedata[seq]['psms'][pool].add(psm_id)
    except KeyError:
        try:
            peptidedata[seq]['psms'][pool] = set()
        except KeyError:
            peptidedata[seq] = {'psms': {pool: set()},
                                'proteins': set()}
        peptidedata[seq]['psms'][pool].add(psm_id)
    if not genecentric:
        protein = (p_acc, cov, desc, gene, assoc_id, len(pgcontent))
    elif genecentric == 'plain':
        protein = (p_acc, cov, desc, gene, assoc_id)
    else:
        protein = (None, cov, desc, gene, assoc_id)
    peptidedata[seq]['proteins'].add(protein)
