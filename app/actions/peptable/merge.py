from app.dataformats import peptable as peptabledata
from app.actions.mergetable import (simple_val_fetch, fill_mergefeature,
                                    parse_NA)
from app.actions.mergetable import create_featuredata_map


def build_peptidetable(pqdb, header, headerfields, isobaric=False,
                       precursor=False, fdr=False, pep=False,
                       peptidedata=False):
    """Fetches peptides and quants from joined lookup table, loops through
    them and when all of a peptides quants/data have been collected, yields
    peptide quant information."""
    peptidedatamap = create_featuredata_map(pqdb, add_record_to_peptidedata, 
                                            get_uniques=False)
    count_psms(peptidedatamap)
    empty_return = lambda x, y, z: {}
    iso_fun = {True: get_isobaric_quant, False: empty_return}[isobaric]
    ms1_fun = {True: get_precursor_quant, False: empty_return}[precursor]
    fdr_fun = {True: get_pep_fdr,
               False: empty_return}[fdr]
    pep_fun = {True: get_peptide_pep,
               False: empty_return}[pep]
    pdata_fun = {True: get_peptide_data, False: empty_return}[peptidedata]
    peptide_sql, sqlfieldmap = pqdb.prepare_mergetable_sql(precursor, isobaric,
                                                           fdr, pep)
    peptides = pqdb.get_merged_features(peptide_sql)
    peptide = next(peptides)
    outpeptide = {peptabledata.HEADER_PEPTIDE: peptide[sqlfieldmap['p_seq']]}
    fill_mergefeature(outpeptide, iso_fun, ms1_fun, empty_return, fdr_fun,
                      pep_fun, pdata_fun, peptide, sqlfieldmap, headerfields,
                      peptidedatamap)
    for peptide in peptides:
        p_seq = peptide[sqlfieldmap['p_seq']]
        if p_seq != outpeptide[peptabledata.HEADER_PEPTIDE]:
            yield parse_NA(outpeptide, header)
            outpeptide = {peptabledata.HEADER_PEPTIDE: p_seq}
        fill_mergefeature(outpeptide, iso_fun, ms1_fun, empty_return, fdr_fun,
                          pep_fun, pdata_fun, peptide, sqlfieldmap,
                          headerfields, peptidedatamap)
    yield parse_NA(outpeptide, header)


def get_isobaric_quant(peptide, sqlmap, headerfields):
    chan = peptide[sqlmap['channel']]
    pool = peptide[sqlmap['isoq_poolname']]
    quant = peptide[sqlmap['isoq_val']]
    return {headerfields['isoquant'][chan][pool]: quant}


def get_precursor_quant(peptide, sqlmap, headerfields):
    return simple_val_fetch(peptide, sqlmap,
                            headerfields['precursorquant'][
                                peptabledata.HEADER_AREA],
                            'preq_poolname', 'preq_val')


def get_pep_fdr(peptide, sqlmap, headerfields):
    return simple_val_fetch(peptide, sqlmap,
                            headerfields['peptidefdr'][
                                peptabledata.HEADER_QVAL],
                            'fdr_poolname', 'fdr_val')


def get_peptide_pep(peptide, sqlmap, headerfields):
    return simple_val_fetch(peptide, sqlmap,
                            headerfields['peptidepep'][
                                peptabledata.HEADER_PEP],
                            'pep_poolname', 'pep_val')


def get_peptide_data(peptide, pdata, headerfields):
    seq = peptide[peptabledata.HEADER_PEPTIDE]
    outdict = {}
    # first the header fields that are pool-dependent
    hfields = [peptabledata.HEADER_NO_PSM,
               ]
    for pool, psms in pdata[seq]['psms'].items():
        pool_values = [psms]
        outdict.update({headerfields['peptidedata'][hfield][pool]: val
                        for (hfield, val) in zip(hfields, pool_values)})
    # Now the pool-independent fields
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
