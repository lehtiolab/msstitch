from app.dataformats import prottable as prottabledata
from app.actions.prottable import info as pinfo
from app.actions.mergetable import (simple_val_fetch, fill_mergefeature,
                                    parse_NA, create_featuredata_map)


def build_proteintable(pqdb, header, headerfields, isobaric=False,
                       precursor=False, probability=False, fdr=False,
                       pep=False, proteindata=False):
    """Fetches proteins and quants from joined lookup table, loops through
    them and when all of a protein's quants have been collected, yields the
    protein quant information."""
    pdmap = create_featuredata_map(pqdb, pinfo.count_peps_psms,
                                   pinfo.add_record_to_proteindata)
    empty_return = lambda x, y, z: {}
    iso_fun = {True: get_isobaric_quant, False: empty_return}[isobaric]
    ms1_fun = {True: get_precursor_quant, False: empty_return}[precursor]
    prob_fun = {True: get_prot_probability,
                False: empty_return}[probability]
    fdr_fun = {True: get_prot_fdr,
               False: empty_return}[fdr]
    pep_fun = {True: get_prot_pep,
               False: empty_return}[pep]
    pdata_fun = {True: get_protein_data, False: empty_return}[proteindata]
    protein_sql, sqlfieldmap = pqdb.prepare_mergetable_sql(precursor, isobaric,
                                                           probability, fdr,
                                                           pep)
    proteins = pqdb.get_merged_features(protein_sql)
    protein = next(proteins)
    outprotein = {prottabledata.HEADER_PROTEIN: protein[sqlfieldmap['p_acc']]}
    fill_mergefeature(outprotein, iso_fun, ms1_fun, prob_fun, fdr_fun, pep_fun,
                      pdata_fun, protein, sqlfieldmap, headerfields,
                      pdmap)
    for protein in proteins:
        p_acc = protein[sqlfieldmap['p_acc']]
        if p_acc != outprotein[prottabledata.HEADER_PROTEIN]:
            yield parse_NA(outprotein, header)
            outprotein = {prottabledata.HEADER_PROTEIN: p_acc}
        fill_mergefeature(outprotein, iso_fun, ms1_fun, prob_fun, fdr_fun,
                          pep_fun, pdata_fun, protein, sqlfieldmap,
                          headerfields, pdmap)
    yield parse_NA(outprotein, header)


def get_protein_data(outprotein, proteindata_map, headerfields):
    p_acc = outprotein[prottabledata.HEADER_PROTEIN]
    return pinfo.get_protein_data(proteindata_map,
                                  p_acc, headerfields)


def get_isobaric_quant(protein, sqlmap, headerfields):
    chan = protein[sqlmap['channel']]
    pool = protein[sqlmap['isoq_poolname']]
    psmfield = protein[sqlmap['isoq_psmsfield']]
    quant = protein[sqlmap['isoq_val']]
    nopsms = protein[sqlmap['isoq_psms']]
    return {headerfields['isoquant'][chan][pool]: quant,
            headerfields['isoquant'][psmfield][pool]: nopsms}


def get_precursor_quant(protein, sqlmap, headerfields):
    return simple_val_fetch(protein, sqlmap,
                            headerfields['precursorquant'][
                                prottabledata.HEADER_AREA],
                            'preq_poolname', 'preq_val')


def get_prot_probability(protein, sqlmap, headerfields):
    return simple_val_fetch(protein, sqlmap,
                            headerfields['probability'][
                                prottabledata.HEADER_PROBABILITY],
                            'prob_poolname', 'prob_val')


def get_prot_fdr(protein, sqlmap, headerfields):
    return simple_val_fetch(protein, sqlmap,
                            headerfields['proteinfdr'][
                                prottabledata.HEADER_QVAL],
                            'fdr_poolname', 'fdr_val')


def get_prot_pep(protein, sqlmap, headerfields):
    return simple_val_fetch(protein, sqlmap,
                            headerfields['proteinpep'][
                                prottabledata.HEADER_PEP],
                            'pep_poolname', 'pep_val')
