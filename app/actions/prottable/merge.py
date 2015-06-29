from app.dataformats import prottable as prottabledata


def build_proteintable(pqdb, header, headerfields, isobaric=False, precursor=False, probability=False, proteindata=False):
    """Fetches proteins and quants from joined lookup table, loops through
    them and when all of a protein's quants have been collected, yields the
    protein quant information."""
    print(isobaric, precursor, probability, proteindata)
    iso_fun = {True: get_isobaric_quant, False: lambda x, y, z: {}}[isobaric]
    ms1_fun = {True: get_precursor_quant, False: lambda x, y, z: {}}[precursor]
    prob_fun = {True: get_prot_probability, False: lambda x, y, z: {}}[probability]
    #pdata_fun = {True: get_proteindata, False: lambda x: {}}[proteindata]
    protein_sql, sqlfieldmap = pqdb.prepare_mergetable_sql(precursor, isobaric, probability)
    proteins = pqdb.get_merged_proteins(protein_sql)
    protein = next(proteins)
    print(protein)
    outprotein = {prottabledata.HEADER_PROTEIN: protein[sqlfieldmap['p_acc']]}
    fill_outprotein(outprotein, iso_fun, ms1_fun, prob_fun, protein, sqlfieldmap, headerfields)
    for protein in proteins:
        if protein[sqlfieldmap['p_acc']] != outprotein[prottabledata.HEADER_PROTEIN]:
            #yield parse_NA(next(add_protein_data([outprotein], pqdb)), header)
            yield parse_NA(outprotein, header)
            outprotein = {prottabledata.HEADER_PROTEIN: protein[sqlfieldmap['p_acc']]}
        fill_outprotein(outprotein, iso_fun, ms1_fun, prob_fun, protein, sqlfieldmap, headerfields)
    #yield parse_NA(next(add_protein_data([outprotein], pqdb)), header)
    yield parse_NA(outprotein, header)


def fill_outprotein(outprotein, iso_fun, ms1_fun, prob_fun, protein, sqlfieldmap, headerfields):
    outprotein.update(iso_fun(protein, sqlfieldmap, headerfields))
    outprotein.update(ms1_fun(protein, sqlfieldmap, headerfields))
    outprotein.update(prob_fun(protein, sqlfieldmap, headerfields))


def add_protein_data(proteins, pgdb):
    """Loops proteins and calls a parsing method to get information
    from a lookup db. Yields proteins with output data"""
    for protein in proteins:
        outprotein = {k: v for k, v in protein.items()}
        protein_acc = protein[prottabledata.HEADER_PROTEIN]
        outprotein.update(get_protein_data(protein_acc, pgdb))
        outprotein = {k: str(v) for k, v in outprotein.items()}
        yield outprotein


#def get_protein_data(protein_acc, pgdb):
#    """Parses protein data that is fetched from the database."""
#    #protein data is ((psm_id, psmseq, fakemaster, all_group_proteins_acc,
#    #                   coverage, description),)
#    protein_data = pgdb.get_protein_data(protein_acc)
#    description = protein_data[0][5]
#    coverage = protein_data[0][4]
#    psmcount = len(set([x[0] for x in protein_data]))
#    pepcount = len(set([x[1] for x in protein_data]))
#    proteincount = len(set([x[3] for x in protein_data]))
#    peptides_master_map = {}
#    for psm in protein_data:
#        try:
#            peptides_master_map[psm[1]].add(psm[2])
#        except KeyError:
#            peptides_master_map[psm[1]] = {psm[2]}
#    unipepcount = len([x for x in peptides_master_map
#                       if len(peptides_master_map[x]) == 1])
#    return {prottabledata.HEADER_DESCRIPTION: description,
#            prottabledata.HEADER_COVERAGE: coverage,
#            prottabledata.HEADER_NO_PROTEIN: proteincount,
#            prottabledata.HEADER_NO_UNIPEP: unipepcount,
#            prottabledata.HEADER_NO_PEPTIDE: pepcount,
#            prottabledata.HEADER_NO_PSM: psmcount,
#            #prottabledata.HEADER_AREA: area,
#            #prottabledata.HEADER_NO_QUANT_PSM: quantcount,
#            #prottabledata.HEADER_CV_QUANT_PSM: quantcv,
#            }
#

def get_isobaric_quant(protein, sqlmap, headerfields):
    chan = protein[sqlmap['channel']]
    pool = protein[sqlmap['isoq_poolname']]
    psmfield = protein[sqlmap['isoq_psmsfield']]
    quant = protein[sqlmap['isoq_val']]
    nopsms = protein[sqlmap['isoq_psms']]
    return {headerfields['isoquant'][chan][pool]: quant,
            headerfields['isoquant'][psmfield][pool]: nopsms}


def get_precursor_quant(protein, sqlmap, headerfields):
    pool = protein[sqlmap['preq_poolname']]
    quantheadfield = headerfields['precursorquant'][prottabledata.HEADER_AREA][pool]
    return {quantheadfield: protein[sqlmap['preq_val']]}


def get_prot_probability(protein, sqlmap, headerfields):
    pool = protein[sqlmap['prob_poolname']]
    headfield = headerfields['probability'][prottabledata.HEADER_PROBABILITY][pool]
    return {headfield: protein[sqlmap['prob_val']]}


def parse_NA(protein, header):
    for field in header:
        try:
            protein[field] = str(protein[field])
        except KeyError:
            protein[field] = 'NA'
    return protein
