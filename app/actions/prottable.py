from app.readers import tsv as tsvreader
from app.dataformats import prottable as prottabledata
from app.dataformats import mzidtsv as mzidtsvdata


def add_ms1_quant_from_top3_mzidtsv(proteins, psms):
    """Collects PSMs with the highes precursor quant values,
    adds sum of the top 3 of these to a protein table"""
    top_ms1_psms = {}
    for psm in psms:
        protacc = psm[mzidtsvdata.HEADER_MASTER_PROT]
        precursor_amount = psm[mzidtsvdata.HEADER_PRECURSOR_QUANT]
        if ';' in protacc or precursor_amount == 'NA':
            continue
        precursor_amount = float(precursor_amount)
        psm_seq = psm[mzidtsvdata.HEADER_PEPTIDE]
        try:
            min_precursor_amount = min(top_ms1_psms[protacc])
        except KeyError:
            top_ms1_psms[protacc] = {-1: None, -2: None,
                                         precursor_amount: psm_seq}
            continue
        else:
            if precursor_amount > min_precursor_amount:
                top_ms1_psms[protacc][precursor_amount] = psm_seq 
                top_ms1_psms[protacc].pop(min_precursor_amount)
    for protein in proteins:
        outprotein = {k: v for k, v in protein.items()}
        try:
            amounts = top_ms1_psms[protein[prottabledata.HEADER_PROTEIN]]
        except KeyError:
            prec_area = 'NA'
        else:
            amounts = [x for x in amounts if x > 0]
            prec_area = sum(amounts) / len(amounts)
        outprotein[prottabledata.HEADER_AREA] = str(prec_area)
        yield outprotein


def get_quantchannels(pqdb):
    quantheader = []
    for fn, chan_name, amnt_psms_name in pqdb.get_quantchannel_headerfields():
        quantheader.append(build_quantchan_header_field(fn, chan_name))
        quantheader.append(build_quantchan_header_field(fn, amnt_psms_name))
    return sorted(quantheader)


def get_precursorquant_headerfields(pqdb):
    return pqdb.get_precursorquant_headerfields()


def build_quantchan_header_field(fn, channame):
    return '{}_{}'.format(fn, channame)


def get_header(oldheader=None, quant_psm_channels=None, addprotein_data=False,
               addprecursor_area=False):
    if oldheader is None:
        header = [prottabledata.HEADER_PROTEIN]
        header.extend(quant_psm_channels)
    else:
        header = oldheader[:]
    if addprotein_data:
        header = get_header_with_proteindata(header)
    if addprecursor_area:
        header = get_header_with_precursorarea(header)
    return header


def get_header_with_proteindata(header):
    ix = header.index(prottabledata.HEADER_PROTEIN) + 1
    new_data = [prottabledata.HEADER_DESCRIPTION,
                prottabledata.HEADER_COVERAGE,
                prottabledata.HEADER_NO_PROTEIN,
                prottabledata.HEADER_NO_UNIPEP,
                prottabledata.HEADER_NO_PEPTIDE,
                prottabledata.HEADER_NO_PSM,
                #prottabledata.HEADER_NO_QUANT_PSM,
                #prottabledata.HEADER_CV_QUANT_PSM,
                ]
    return header[:ix] + new_data + header[ix:]


def get_header_with_precursorarea(header):
    ix = header.index(prottabledata.HEADER_PROTEIN) + 1
    return header[:ix] + [prottabledata.HEADER_AREA] + header[ix:]


def build_quanted_proteintable(pqdb, header):
    """Fetches proteins and quants from joined lookup table, loops through
    them and when all of a protein's quants have been collected, yields the
    protein quant information."""
    proteins = pqdb.get_quanted_proteins()
    protein = next(proteins)
    outprotein = {prottabledata.HEADER_PROTEIN: protein[0],
                  protein[1]: protein[2]}
    for protein in proteins:
        if protein[0] != outprotein[prottabledata.HEADER_PROTEIN]:
            yield parse_NA(next(add_protein_data([outprotein], pqdb)), header)
            outprotein = {prottabledata.HEADER_PROTEIN: protein[0]}
        quantheadfield = build_quantchan_header_field(protein[2], protein[1])
        amntpsm_headfld = build_quantchan_header_field(protein[2], protein[3])
        outprotein[quantheadfield] = protein[4]
        outprotein[amntpsm_headfld] = protein[5]
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


def parse_NA(protein, header):
    for field in header:
        try:
            protein[field]
        except KeyError:
            protein[field] = 'NA'
    return protein
