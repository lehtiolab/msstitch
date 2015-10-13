from app.dataformats import prottable as prottabledata
from app.actions.prottable import info as pinfo
from app.actions.mergetable import (simple_val_fetch, fill_mergefeature,
                                    parse_NA, create_featuredata_map)


def build_proteintable(pqdb, header, headerfields, isobaric=False,
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
    psms_fun = empty_return
    pdata_fun = {True: get_protein_data_genecentric,
                 False: get_protein_data}[genecentric]
    protein_sql, sqlfieldmap = pqdb.prepare_mergetable_sql(precursor, isobaric,
                                                           probability, fdr,
                                                           pep)
    proteins = pqdb.get_merged_features(protein_sql)
    protein = next(proteins)
    outprotein = {prottabledata.HEADER_PROTEIN: protein[sqlfieldmap['p_acc']]}
    check_prot = {k: v for k, v in outprotein.items()}
    fill_mergefeature(outprotein, iso_fun, ms1_fun, prob_fun, fdr_fun, pep_fun,
                      psms_fun, pdata_fun, protein, sqlfieldmap, headerfields,
                      pdmap)
    for protein in proteins:
        p_acc = protein[sqlfieldmap['p_acc']]
        if p_acc != outprotein[prottabledata.HEADER_PROTEIN]:
            if outprotein != check_prot:
                yield parse_NA(outprotein, header)
            outprotein = {prottabledata.HEADER_PROTEIN: p_acc}
            check_prot = {k: v for k, v in outprotein.items()}
        fill_mergefeature(outprotein, iso_fun, ms1_fun, prob_fun, fdr_fun,
                          pep_fun, psms_fun, pdata_fun, protein, sqlfieldmap,
                          headerfields, pdmap)
    yield parse_NA(outprotein, header)


def get_protein_data_genecentric(outprotein, proteindata_map, headerfields):
    p_acc = outprotein[prottabledata.HEADER_PROTEIN]
    return pinfo.get_protein_data_genecentric(proteindata_map, p_acc,
                                              headerfields)


def get_protein_data(outprotein, proteindata_map, headerfields):
    p_acc = outprotein[prottabledata.HEADER_PROTEIN]
    return pinfo.get_protein_data_pgrouped(proteindata_map, p_acc,
                                           headerfields)


def get_isobaric_quant(protein, sqlmap, headerfields):
    chan = protein[sqlmap['channel']]
    pool = protein[sqlmap['set_name']]
    psmfield = protein[sqlmap['isoq_psmsfield']]
    quant = protein[sqlmap['isoq_val']]
    nopsms = protein[sqlmap['isoq_psms']]
    if quant is None:
        return {}
    return {headerfields['isoquant'][chan][pool]: quant,
            headerfields['isoquant'][psmfield][pool]: nopsms}


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
