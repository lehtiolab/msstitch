from app.dataformats import prottable as prottabledata
from app.actions.prottable import info as pdatagenerator


def build_proteintable(pqdb, header, headerfields, isobaric=False,
                       precursor=False, probability=False, proteindata=False):
    """Fetches proteins and quants from joined lookup table, loops through
    them and when all of a protein's quants have been collected, yields the
    protein quant information."""
    proteindatamap = pdatagenerator.create_proteindata_map(pqdb)
    iso_fun = {True: get_isobaric_quant, False: lambda x, y, z: {}}[isobaric]
    ms1_fun = {True: get_precursor_quant, False: lambda x, y, z: {}}[precursor]
    prob_fun = {True: get_prot_probability,
                False: lambda x, y, z: {}}[probability]
    pdata_fun = {True: get_protein_data, False: lambda x: {}}[proteindata]
    protein_sql, sqlfieldmap = pqdb.prepare_mergetable_sql(precursor, isobaric,
                                                           probability)
    proteins = pqdb.get_merged_proteins(protein_sql)
    protein = next(proteins)
    outprotein = {prottabledata.HEADER_PROTEIN: protein[sqlfieldmap['p_acc']]}
    fill_outprotein(outprotein, iso_fun, ms1_fun, prob_fun, pdata_fun, protein,
                    sqlfieldmap, headerfields, proteindatamap)
    for protein in proteins:
        p_acc = protein[sqlfieldmap['p_acc']]
        if p_acc != outprotein[prottabledata.HEADER_PROTEIN]:
            yield parse_NA(outprotein, header)
            outprotein = {prottabledata.HEADER_PROTEIN: p_acc}
        fill_outprotein(outprotein, iso_fun, ms1_fun, prob_fun, pdata_fun,
                        protein, sqlfieldmap, headerfields, proteindatamap)
    yield parse_NA(outprotein, header)


def fill_outprotein(outprotein, iso_fun, ms1_fun, prob_fun, pdata_fun, protein,
                    sqlfieldmap, headerfields, proteindata_map):
    outprotein.update(iso_fun(protein, sqlfieldmap, headerfields))
    outprotein.update(ms1_fun(protein, sqlfieldmap, headerfields))
    outprotein.update(prob_fun(protein, sqlfieldmap, headerfields))
    outprotein.update(pdata_fun(outprotein, proteindata_map, headerfields))


def get_protein_data(outprotein, proteindata_map, headerfields):
    p_acc = outprotein[prottabledata.HEADER_PROTEIN]
    return pdatagenerator.get_protein_data(proteindata_map,
                                           p_acc, headerfields)


def get_isobaric_quant(protein, sqlmap, headerfields):
    chan = protein[sqlmap['channel']]
    pool = protein[sqlmap['isoq_poolname']]
    psmfield = protein[sqlmap['isoq_psmsfield']]
    quant = protein[sqlmap['isoq_val']]
    nopsms = protein[sqlmap['isoq_psms']]
    return {headerfields['isoquant'][chan][pool]: quant,
            headerfields['isoquant'][psmfield][pool]: nopsms}


def simple_val_fetch(protein, sqlmap, headerfields, poolkey, fieldkey, valkey):
    pool = protein[sqlmap[poolkey]]
    hfield = headerfields[fieldkey][prottabledata.HEADER_AREA][pool]
    return {hfield: protein[sqlmap[valkey]]}


def get_precursor_quant(protein, sqlmap, headerfields):
    return simple_val_fetch(protein, sqlmap, headerfields, 'preq_poolname',
                            'precursorquant', 'preq_val')


def get_prot_probability(protein, sqlmap, headerfields):
    return simple_val_fetch(protein, sqlmap, headerfields, 'prob_poolname',
                            'probability', 'prob_val')


def parse_NA(protein, header):
    for field in header:
        try:
            protein[field] = str(protein[field])
        except KeyError:
            protein[field] = 'NA'
    return protein
