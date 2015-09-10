from app.dataformats import peptable as peptabledata
from app.actions.mergetable import (simple_val_fetch, fill_mergefeature,
                                    parse_NA)
#from app.actions.prottable import info as pdatagenerator


def build_peptable(pqdb, header, headerfields, isobaric=False, precursor=False,
                   fdr=False, pep=False, peptidedata=False):
    """Fetches peptides and quants from joined lookup table, loops through
    them and when all of a peptides quants/data have been collected, yields
    peptide quant information."""
    #peptidedatamap = pdatagenerator.create_peptidedata_map(pqdb)
    peptidedatamap = None
    empty_return = lambda x, y, z: {}
    iso_fun = {True: get_isobaric_quant, False: empty_return}[isobaric]
    ms1_fun = {True: get_precursor_quant, False: empty_return}[precursor]
    fdr_fun = {True: get_pep_fdr,
               False: empty_return}[fdr]
    pep_fun = {True: get_peptide_pep,
               False: empty_return}[pep]
    pdata_fun = empty_return
    # {True: get_peptide_data, False: empty_return}[peptidedata]
    peptide_sql, sqlfieldmap = pqdb.prepare_mergetable_sql(precursor, isobaric,
                                                           fdr, pep)
    peptides = pqdb.get_merged_peptides(peptide_sql)
    peptide = next(peptides)
    outpeptide = {peptabledata.HEADER_PEPTIDE: peptide[sqlfieldmap['p_acc']]}
    fill_mergefeature(outpeptide, iso_fun, ms1_fun, fdr_fun, pep_fun,
                      pdata_fun, peptide, sqlfieldmap, headerfields,
                      peptidedatamap)
    for peptide in peptides:
        p_acc = peptide[sqlfieldmap['p_acc']]
        if p_acc != outpeptide[peptabledata.HEADER_PROTEIN]:
            yield parse_NA(outpeptide, header)
            outpeptide = {peptabledata.HEADER_PROTEIN: p_acc}
        fill_mergefeature(outpeptide, iso_fun, ms1_fun, fdr_fun,
                          pep_fun, pdata_fun, peptide, sqlfieldmap,
                          headerfields, peptidedatamap)
    yield parse_NA(outpeptide, header)


def get_isobaric_quant(protein, sqlmap, headerfields):
    chan = protein[sqlmap['channel']]
    pool = protein[sqlmap['isoq_poolname']]
    quant = protein[sqlmap['isoq_val']]
    return {headerfields['isoquant'][chan][pool]: quant}


def get_precursor_quant(peptide, sqlmap, headerfields):
    return simple_val_fetch(peptide, sqlmap,
                            headerfields['precursorquant'][peptabledata.HEADER_AREA],
                            'preq_poolname', 'preq_val')


def get_pep_fdr(peptide, sqlmap, headerfields):
    return simple_val_fetch(peptide, sqlmap,
                            headerfields['peptidefdr'][peptabledata.HEADER_QVAL],
                            'fdr_poolname', 'fdr_val')


def get_peptide_pep(peptide, sqlmap, headerfields):
    return simple_val_fetch(peptide, sqlmap,
                            headerfields['peptidepep'][peptabledata.HEADER_PEP],
                            'pep_poolname', 'pep_val')
