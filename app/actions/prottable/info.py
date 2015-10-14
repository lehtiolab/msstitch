from app.dataformats import prottable as prottabledata
from app.actions.mergetable import create_featuredata_map


def add_record_to_proteindata(proteindata, p_acc, pool, psmdata, genecentric,
                              pgcontentmap=None):
    """Fill function for create_featuredata_map"""
    seq, psm_id = psmdata[2], psmdata[3]
    if not genecentric:
        desc, cov = psmdata[4], psmdata[5]
        pgcontent = pgcontentmap[p_acc]
    else:
        cov = None
    try:
        proteindata[p_acc]['pools'][pool]['psms'].add(psm_id)
    except KeyError:
        emptyinfo = {'psms': set(), 'peptides': set(), 'unipeps': 0}
        try:
            proteindata[p_acc]['pools'][pool] = emptyinfo
        except KeyError:
            proteindata[p_acc] = {'pools': {pool: emptyinfo}}
            if cov is not None:
                proteindata[p_acc].update({'desc': desc, 'cov': cov, 
                                           'proteins': pgcontent})
    proteindata[p_acc]['pools'][pool]['psms'].add(psm_id)
    proteindata[p_acc]['pools'][pool]['peptides'].add(seq)


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
                                         pool_to_output=pool_to_output)
    dataget_fun = {True: get_protein_data_genecentric,
                   False: get_protein_data_pgrouped}[genecentric]
    for protein in proteins:
        outprotein = {k: v for k, v in protein.items()}
        protein_acc = protein[prottabledata.HEADER_PROTEIN]
        if not protein_acc in proteindata:
            continue
        outprotein.update(dataget_fun(proteindata, protein_acc, headerfields))
        outprotein = {k: str(v) for k, v in outprotein.items()}
        yield outprotein


def get_protein_data_genecentric(proteindata, p_acc, headerfields):
    return get_protein_data_base(proteindata, p_acc, headerfields)


def get_protein_data_pgrouped(proteindata, p_acc, headerfields):
    """Parses protein data for a certain protein into tsv output
    dictionary"""
    report = get_protein_data_base(proteindata, p_acc, headerfields)
    return get_cov_descriptions(proteindata, p_acc, report)


def get_protein_data_base(proteindata, p_acc, headerfields):
    proteincount = 'na'
    outdict = {}
    hfields = [prottabledata.HEADER_NO_UNIPEP,
               prottabledata.HEADER_NO_PEPTIDE,
               prottabledata.HEADER_NO_PSM,
               ]
    for pool, pdata in proteindata[p_acc]['pools'].items():
        pool_values = [pdata['unipeps'], pdata['peptides'], pdata['psms']]
        outdict.update({headerfields['proteindata'][hfield][pool]: val
                        for (hfield, val) in zip(hfields, pool_values)})
    outdict[prottabledata.HEADER_NO_PROTEIN] = proteincount
    return outdict


def get_cov_descriptions(proteindata, p_acc, report):
    try:
        report.update({prottabledata.HEADER_DESCRIPTION:
                       proteindata[p_acc]['desc'],
                       prottabledata.HEADER_COVERAGE:
                       proteindata[p_acc]['cov'],
                       prottabledata.HEADER_CONTENTPROT:
                       ','.join(proteindata[p_acc]['proteins']),
                       prottabledata.HEADER_NO_PROTEIN:
                       len(proteindata[p_acc]['proteins'])
                       })
    except KeyError:
        # In case database is built without a FASTA file there is no coverage
        # info available.
        pass
    return report
