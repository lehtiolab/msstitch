from app.dataformats import prottable as prottabledata


def build_proteintable(pqdb, header, isobaric=False, precursor=False, probability=False):
    """Fetches proteins and quants from joined lookup table, loops through
    them and when all of a protein's quants have been collected, yields the
    protein quant information."""
    iso_quant_map = {True: get_isobaric_quant, False: lambda x, y, z: {}}
    ms1_quant_map = {True: get_precursor_quant, False: lambda x, y, z: {}}
    prob_map = {True: get_prot_probability, False: lambda x, y, z: {}}
    inv_prottable_map = get_inverted_prottable_map(pqdb)
    protein_sql, sqlfieldmap = pqdb.prepare_mergetable_sql(precursor, isobaric, probability)
    proteins = pqdb.get_merged_proteins(protein_sql)
    protein = next(proteins)
    outprotein = {prottabledata.HEADER_PROTEIN: protein[sqlfieldmap['p_acc']]}
    outprotein.update(iso_quant_map[isobaric](protein, sqlfieldmap, inv_prottable_map))
    outprotein.update(ms1_quant_map[precursor](protein, sqlfieldmap, inv_prottable_map))
    for protein in proteins:
        if protein[sqlfieldmap['p_acc']] != outprotein[prottabledata.HEADER_PROTEIN]:
            yield parse_NA(next(add_protein_data([outprotein], pqdb)), header)
            outprotein = {prottabledata.HEADER_PROTEIN: protein[sqlfieldmap['p_acc']]}
        outprotein.update(iso_quant_map[isobaric](protein, sqlfieldmap, inv_prottable_map))
        outprotein.update(ms1_quant_map[precursor](protein, sqlfieldmap, inv_prottable_map))
    yield parse_NA(next(add_protein_data([outprotein], pqdb)), header)


def add_protein_data(proteins, pgdb):
    """Loops proteins and calls a parsing method to get information
    from a lookup db. Yields proteins with output data"""
    for protein in proteins:
        outprotein = {k: v for k, v in protein.items()}
        protein_acc = protein[prottabledata.HEADER_PROTEIN]
        outprotein.update(get_protein_data(protein_acc, pgdb))
        outprotein = {k: str(v) for k, v in outprotein.items()}
        yield outprotein


def get_protein_data(protein_acc, pgdb):
    """Parses protein data that is fetched from the database."""
    #protein data is ((psm_id, psmseq, fakemaster, all_group_proteins_acc,
    #                   coverage, description),)
    protein_data = pgdb.get_protein_data(protein_acc)
    description = protein_data[0][5]
    coverage = protein_data[0][4]
    psmcount = len(set([x[0] for x in protein_data]))
    pepcount = len(set([x[1] for x in protein_data]))
    proteincount = len(set([x[3] for x in protein_data]))
    peptides_master_map = {}
    for psm in protein_data:
        try:
            peptides_master_map[psm[1]].add(psm[2])
        except KeyError:
            peptides_master_map[psm[1]] = {psm[2]}
    unipepcount = len([x for x in peptides_master_map
                       if len(peptides_master_map[x]) == 1])
    return {prottabledata.HEADER_DESCRIPTION: description,
            prottabledata.HEADER_COVERAGE: coverage,
            prottabledata.HEADER_NO_PROTEIN: proteincount,
            prottabledata.HEADER_NO_UNIPEP: unipepcount,
            prottabledata.HEADER_NO_PEPTIDE: pepcount,
            prottabledata.HEADER_NO_PSM: psmcount,
            #prottabledata.HEADER_AREA: area,
            #prottabledata.HEADER_NO_QUANT_PSM: quantcount,
            #prottabledata.HEADER_CV_QUANT_PSM: quantcv,
            }


def get_isobaric_quant(protein, sqlmap, prottable_map):
    chan = protein[sqlmap['channel']]
    fn = prottable_map[protein[sqlmap['isoq_fnid']]]
    psmfield = protein[sqlmap['isoq_psmsfield']]
    quant = protein[sqlmap['isoq_val']]
    nopsms = protein[sqlmap['isoq_psms']]
    quantheadfield = build_quantchan_header_field(fn, chan)
    amntpsm_headfld = build_quantchan_header_field(fn, psmfield)
    return {quantheadfield: quant, amntpsm_headfld: nopsms}


def get_precursor_quant(protein, sqlmap, prottable_map):
    fn = prottable_map[protein[sqlmap['preq_fnid']]]
    quantheadfield = build_quantchan_header_field(fn, prottabledata.HEADER_AREA)
    return {quantheadfield: protein[sqlmap['preq_val']]}


def get_precursorquant_headerfields(pqdb):
    prottable_map = get_inverted_prottable_map(pqdb)
    fnids = pqdb.get_precursorquant_headerfields()
    return [prottable_map[fnid[0]] for fnid in fnids]


def build_quantchan_header_field(fn, channame):
    # FIXME should be in header module
    return '{}_{}'.format(fn, channame)


def get_inverted_prottable_map(pqdb):
    prottable_map = pqdb.get_protein_table_map()
    return {v: k for k, v in prottable_map.items()}

def get_quantchannels(pqdb):
    # FIXME should be in header module
    quantheader = []
    prottable_map = get_inverted_prottable_map(pqdb)
    for fnid, chan_name, amnt_psms_name in pqdb.get_quantchannel_headerfields():
        quantheader.append(build_quantchan_header_field(prottable_map[fnid], chan_name))
        quantheader.append(build_quantchan_header_field(prottable_map[fnid], amnt_psms_name))
    return sorted(quantheader)


def parse_NA(protein, header):
    for field in header:
        try:
            protein[field]
        except KeyError:
            protein[field] = 'NA'
    return protein
