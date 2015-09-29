from app.dataformats import prottable as prottabledata
from app.actions.mergetable import create_featuredata_map


def add_record_to_proteindata(proteindata, p_acc, pool, psmdata):
    seq, psm_id = psmdata[2], psmdata[3]
    desc, cov = psmdata[4], psmdata[5]
    try:
        proteindata[p_acc]['pools'][pool]['psms'].add(psm_id)
    except KeyError:
        emptyinfo = {'psms': set(), 'peptides': set(), 'proteins': set(),
                     'unipeps': 0}
        try:
            proteindata[p_acc]['pools'][pool] = emptyinfo
        except KeyError:
            proteindata[p_acc] = {'pools': {pool: emptyinfo}}
            if cov is not None:
                proteindata[p_acc].update({'desc': desc, 'cov': cov})
    proteindata[p_acc]['pools'][pool]['psms'].add(psm_id)
    proteindata[p_acc]['pools'][pool]['peptides'].add(seq)
    #proteindata[p_acc][pool]['proteins'].add(pg_content)


def count_peps_psms(proteindata, p_acc, pool):
    data = proteindata[p_acc]['pools'][pool]
    proteindata[p_acc]['pools'][pool]['psms'] = len(data['psms'])
    #proteindata[prot][pool]['proteins'] = len(data['proteins'])
    proteindata[p_acc]['pools'][pool]['peptides'] = len(data['peptides'])


def add_protein_data(proteins, pgdb, headerfields, pool_to_output=False):
    """First creates a map with all master proteins with data,
    then outputs protein data dicts for rows of a tsv. If a pool
    is given then only output for that pool will be shown in the
    protein table."""
    proteindata = create_featuredata_map(pgdb, add_record_to_proteindata,
                                         count_peps_psms, pool_to_output)
    for protein in proteins:
        outprotein = {k: v for k, v in protein.items()}
        protein_acc = protein[prottabledata.HEADER_PROTEIN]
        if not protein_acc in proteindata:
            continue
        outprotein.update(get_protein_data(proteindata, protein_acc,
                                           headerfields))
        outprotein = {k: str(v) for k, v in outprotein.items()}
        yield outprotein


def get_protein_data(proteindata, p_acc, headerfields):
    """Parses protein data for a certain protein into tsv output
    dictionary"""
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
    try:
        outdict.update({prottabledata.HEADER_DESCRIPTION:
                        proteindata[p_acc]['desc'],
                        prottabledata.HEADER_COVERAGE:
                        proteindata[p_acc]['cov'],
                        })
    except KeyError:
        # In case database is built without a FASTA file there is no coverage
        # info available.
        pass
    return outdict
