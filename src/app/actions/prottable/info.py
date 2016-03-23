from app.dataformats import prottable as prottabledata
from app.actions.proteindata import (add_record_to_proteindata,
                                     create_featuredata_map)


def count_peps_psms(proteindata, p_acc, pool):
    data = proteindata[p_acc]['pools'][pool]
    proteindata[p_acc]['pools'][pool]['psms'] = len(data['psms'])
    #proteindata[prot][pool]['proteins'] = len(data['proteins'])
    proteindata[p_acc]['pools'][pool]['peptides'] = len(data['peptides'])


def add_protein_data(proteins, pgdb, headerfields, genecentric=False,
                     pool_to_output=False):
    """First creates a map with all master proteins with data,
    then outputs protein data dicts for rows of a tsv. If a pool
    is given then only output for that pool will be shown in the
    protein table."""
    proteindata = create_featuredata_map(pgdb, genecentric=genecentric,
                                         fill_fun=add_record_to_proteindata,
                                         count_fun=count_peps_psms,
                                         pool_to_output=pool_to_output,
                                         get_uniques=True)
    dataget_fun = {True: get_protein_data_genecentric,
                   False: get_protein_data_pgrouped}[genecentric is not False]
    for protein in proteins:
        outprotein = {k: v for k, v in protein.items()}
        protein_acc = protein[prottabledata.HEADER_PROTEIN]
        outprotein.update(dataget_fun(proteindata, protein_acc, headerfields))
        outprotein = {k: str(v) for k, v in outprotein.items()}
        yield outprotein


def get_protein_data_genecentric(proteindata, p_acc, headerfields):
    return get_protein_data_base(proteindata, p_acc, headerfields)


def get_protein_data_pgrouped(proteindata, p_acc, headerfields):
    """Parses protein data for a certain protein into tsv output
    dictionary"""
    report = get_protein_data_base(proteindata, p_acc, headerfields)
    return get_cov_protnumbers(proteindata, p_acc, report)


def get_headerfieldtext(headerfields, hfieldtype, pool):
    try:
        text = headerfields['proteindata'][hfieldtype][pool]
    except KeyError:
        text = headerfields['proteindata'][hfieldtype][None]
    return text


def get_protein_data_base(proteindata, p_acc, headerfields):
    proteincount = 'na'
    outdict = {}
    hfields = [prottabledata.HEADER_NO_UNIPEP,
               prottabledata.HEADER_NO_PEPTIDE,
               prottabledata.HEADER_NO_PSM,
               ]
    for pool, pdata in proteindata[p_acc]['pools'].items():
        pool_values = [pdata['unipeps'], pdata['peptides'], pdata['psms']]
        outdict.update({get_headerfieldtext(headerfields, hfield, pool): val
                        for (hfield, val) in zip(hfields, pool_values)})
    outdict.update({prottabledata.HEADER_NO_PROTEIN: proteincount,
                    prottabledata.HEADER_DESCRIPTION:
                    proteindata[p_acc]['desc']})
    for field, pdfield in zip([prottabledata.HEADER_GENE,
                               prottabledata.HEADER_ASSOCIATED],
                              ['gene', 'aid']):
        try:
            outdict[field] = ';'.join(proteindata[p_acc][pdfield])
        except TypeError:
            pass
    return outdict


def get_cov_protnumbers(proteindata, p_acc, report):
    try:
        report[prottabledata.HEADER_COVERAGE] = proteindata[p_acc]['cov']
    except KeyError:
        pass
    if 'proteins' in proteindata[p_acc]:
        report.update({prottabledata.HEADER_CONTENTPROT:
                       ','.join(proteindata[p_acc]['proteins']),
                       prottabledata.HEADER_NO_PROTEIN:
                       len(proteindata[p_acc]['proteins'])})
    return report
