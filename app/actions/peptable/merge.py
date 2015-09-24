from app.dataformats import peptable as peptabledata
from app.actions.mergetable import (simple_val_fetch, fill_mergefeature,
                                    parse_NA)
from app.actions.mergetable import create_featuredata_map


def build_peptidetable(pqdb, header, headerfields, isobaric=False,
                       precursor=False, fdr=False, pep=False,
                       nopsms=False, proteindata=False):
    """Fetches peptides and quants from joined lookup table, loops through
    them and when all of a peptides quants/data have been collected, yields
    peptide quant information."""
    peptidedatamap = False
    if nopsms or proteindata:
        peptidedatamap = create_featuredata_map(pqdb, add_record_to_peptidedata,
                                                get_uniques=False)
    if nopsms:
        count_psms(peptidedatamap)
    empty_return = lambda x, y, z: {}
    iso_fun = {True: get_isobaric_quant, False: empty_return}[isobaric]
    ms1_fun = {True: get_precursor_quant, False: empty_return}[precursor]
    fdr_fun = {True: get_pep_fdr,
               False: empty_return}[fdr]
    pep_fun = {True: get_peptide_pep,
               False: empty_return}[pep]
    psms_fun = {True: get_no_psms, False: empty_return}[nopsms]
    pdata_fun = {True: get_protein_data, False: empty_return}[proteindata]
    peptide_sql, sqlfieldmap = pqdb.prepare_mergetable_sql(precursor, isobaric,
                                                           fdr, pep)
    peptides = pqdb.get_merged_features(peptide_sql)
    peptide = next(peptides)
    outpeptide = {peptabledata.HEADER_PEPTIDE: peptide[sqlfieldmap['p_seq']]}
    check_pep = {k: v for k, v in outpeptide.items()}
    fill_mergefeature(outpeptide, iso_fun, ms1_fun, empty_return, fdr_fun,
                      pep_fun, psms_fun, pdata_fun, peptide, sqlfieldmap,
                      headerfields, peptidedatamap)
    for peptide in peptides:
        p_seq = peptide[sqlfieldmap['p_seq']]
        if p_seq != outpeptide[peptabledata.HEADER_PEPTIDE]:
            if outpeptide != check_pep:
                yield parse_NA(outpeptide, header)
            outpeptide = {peptabledata.HEADER_PEPTIDE: p_seq}
            check_pep = {k: v for k, v in outpeptide.items()}
        fill_mergefeature(outpeptide, iso_fun, ms1_fun, empty_return, fdr_fun,
                          pep_fun, psms_fun, pdata_fun, peptide, sqlfieldmap,
                          headerfields, peptidedatamap)
    yield parse_NA(outpeptide, header)


def get_isobaric_quant(peptide, sqlmap, headerfields):
    chan = peptide[sqlmap['channel']]
    pool = peptide[sqlmap['set_name']]
    quant = peptide[sqlmap['isoq_val']]
    if quant is None:
        return {}
    return {headerfields['isoquant'][chan][pool]: quant}


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
    seq = peptide[peptabledata.HEADER_PEPTIDE]
    outdict = {}
    for idx, key in enumerate([peptabledata.HEADER_PROTEINS,
                               peptabledata.HEADER_DESCRIPTIONS,
                               peptabledata.HEADER_COVERAGES]):
        outdict[key] = ';'.join([str(x[idx]) for x in pdata[seq]['proteins']])
    return outdict


def count_psms(pdata):
    for seq in pdata:
        for pool in pdata[seq]['psms']:
            pdata[seq]['psms'][pool] = len(pdata[seq]['psms'][pool])


def add_record_to_peptidedata(peptidedata, p_acc, pool, psmdata):
    seq, psm_id, desc, cov = psmdata[2], psmdata[3], psmdata[4], psmdata[5]
    try:
        peptidedata[seq]['psms'][pool].add(psm_id)
    except KeyError:
        try:
            peptidedata[seq]['psms'][pool] = set()
        except KeyError:
            peptidedata[seq] = {'psms': {pool: set()},
                                'proteins': set()}
        peptidedata[seq]['psms'][pool].add(psm_id)
    peptidedata[seq]['proteins'].add((p_acc, desc, cov))
