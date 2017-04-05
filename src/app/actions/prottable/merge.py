from app.dataformats import prottable as prottabledata
from app.actions.prottable import info as pinfo
from app.actions.mergetable import (simple_val_fetch, fill_mergefeature,
                                    get_isobaric_quant)

from app.actions.proteindata import create_featuredata_map


def build_proteintable(pqdb, headerfields, mergecutoff, isobaric=False,
                       precursor=False, probability=False, fdr=False,
                       pep=False, genecentric=False):
    """Fetches proteins and quants from joined lookup table, loops through
    them and when all of a protein's quants have been collected, yields the
    protein quant information."""
    pdmap = create_featuredata_map(pqdb, genecentric=genecentric,
                                   fill_fun=pinfo.add_record_to_proteindata,
                                   count_fun=pinfo.count_peps_psms,
                                   get_uniques=True)
    empty_return = lambda x, y, z: {}
    iso_fun = {True: get_isobaric_quant, False: empty_return}[isobaric]
    ms1_fun = {True: get_precursor_quant, False: empty_return}[precursor]
    prob_fun = {True: get_prot_probability,
                False: empty_return}[probability]
    fdr_fun = {True: get_prot_fdr,
               False: empty_return}[fdr]
    pep_fun = {True: get_prot_pep,
               False: empty_return}[pep]
    pdata_fun = {True: get_protein_data_genecentric,
                 False: get_protein_data}[genecentric is not False]
    protein_sql, sqlfieldmap = pqdb.prepare_mergetable_sql(precursor, isobaric,
                                                           probability, fdr,
                                                           pep)
    proteins = pqdb.get_merged_features(protein_sql)
    protein = next(proteins)
    outprotein = {prottabledata.HEADER_PROTEIN: protein[sqlfieldmap['p_acc']]}
    check_prot = {k: v for k, v in outprotein.items()}
    if not mergecutoff or protein_pool_fdr_cutoff(protein, sqlfieldmap,
                                                  mergecutoff):
        fill_mergefeature(outprotein, iso_fun, ms1_fun, prob_fun, fdr_fun,
                          pep_fun, pdata_fun, protein, sqlfieldmap,
                          headerfields, pdmap)
    for protein in proteins:
        if mergecutoff and not protein_pool_fdr_cutoff(protein, sqlfieldmap,
                                                       mergecutoff):
            continue
        p_acc = protein[sqlfieldmap['p_acc']]
        if p_acc != outprotein[prottabledata.HEADER_PROTEIN]:
            if outprotein != check_prot:
                # Is this here to avoid outputting empty proteins?
                yield outprotein
            outprotein = {prottabledata.HEADER_PROTEIN: p_acc}
            check_prot = {k: v for k, v in outprotein.items()}
        fill_mergefeature(outprotein, iso_fun, ms1_fun, prob_fun, fdr_fun,
                          pep_fun, pdata_fun, protein, sqlfieldmap,
                          headerfields, pdmap)
    yield outprotein


def protein_pool_fdr_cutoff(sql_entry, sqlmap, fdrcutoff):
    fdr = sql_entry[sqlmap['fdr_val']]
    warning = False
    if fdr is False:
        warning = True
    try:
        fdr = float(fdr)
    except (TypeError, ValueError):
        warning = True
    if warning:
        print('WARNING, filtering merge table on FDR but protein {} in '
              'set {} has FDR value "{}" in '
              'lookup.'.format(sql_entry[sqlmap['p_acc']],
                               sql_entry[sqlmap['set_name']], fdr))
        return False
    return fdr < fdrcutoff


def get_protein_data_genecentric(outprotein, proteindata_map, headerfields):
    p_acc = outprotein[prottabledata.HEADER_PROTEIN]
    return pinfo.get_protein_data_genecentric(proteindata_map, p_acc,
                                              headerfields)


def get_protein_data(outprotein, proteindata_map, headerfields):
    p_acc = outprotein[prottabledata.HEADER_PROTEIN]
    return pinfo.get_protein_data_pgrouped(proteindata_map, p_acc,
                                           headerfields)


def get_precursor_quant(protein, sqlmap, headerfields):
    return simple_val_fetch(protein, sqlmap,
                            headerfields['precursorquant'][
                                prottabledata.HEADER_AREA], 'preq_val')


def get_prot_probability(protein, sqlmap, headerfields):
    return simple_val_fetch(protein, sqlmap,
                            headerfields['probability'][
                                prottabledata.HEADER_PROBABILITY], 'prob_val')


def get_prot_fdr(protein, sqlmap, headerfields):
    return simple_val_fetch(protein, sqlmap,
                            headerfields['proteinfdr'][
                                prottabledata.HEADER_QVAL], 'fdr_val')


def get_prot_pep(protein, sqlmap, headerfields):
    return simple_val_fetch(protein, sqlmap,
                            headerfields['proteinpep'][
                                prottabledata.HEADER_PEP], 'pep_val')
