from app.readers import fasta
from app.dataformats import prottable as prottabledata


def get_header_with_proteindata(header):
    ix = header.index(prottabledata.HEADER_PROTEIN) + 1
    new_data = [prottabledata.HEADER_DESCRIPTION,
                prottabledata.HEADER_COVERAGE,
                prottabledata.HEADER_NO_PROTEIN,
                prottabledata.HEADER_NO_UNIPEP,
                prottabledata.HEADER_NO_PEPTIDE,
                prottabledata.HEADER_NO_PSM,
                #prottabledata.HEADER_AREA,
                #prottabledata.HEADER_NO_QUANT_PSM,
                #prottabledata.HEADER_CV_QUANT_PSM,
                ]
    return header[:ix] + new_data + header[ix:]


def add_protein_data(proteins, header, pgdb):
    for protein in proteins:
        protein_acc = protein[prottabledata.HEADER_PROTEIN]
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


def collect_descriptions(pgdb, fastafn):
    pgdb.add_tables()
    protein_descriptions = fasta.get_proteins_descriptions(fastafn)
    pgdb.store_descriptions(protein_descriptions)
