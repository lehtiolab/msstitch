from app.dataformats import prottable as prottabledata
from app.readers import tsv as reader


def add_isoquant_data(proteins, quantproteins, quantfields):
    """Runs through a protein table and adds quant data from ANOTHER protein
    table that contains that data."""
    quant_proteinmap = get_quantprotein_map(quantproteins, quantfields)
    for protein in proteins:
        prot_acc = protein[prottabledata.HEADER_PROTEIN]
        outprotein = {k: v for k, v in protein.items()}
        try:
            outprotein.update(quant_proteinmap[prot_acc])
        except KeyError:
            outprotein.update({field: 'NA' for field in quantfields})
        print(outprotein)
        yield outprotein


def get_quantprotein_map(proteins, quantfields):
    """Runs through proteins that are in a quanted protein table, extracts
    and maps their information based on the quantfields list input.
    Map is a dict with protein_accessions as keys."""
    qmap = {}
    for protein in proteins:
        prot_acc = protein.pop(prottabledata.HEADER_PROTEIN)
        qmap[prot_acc] = {qf: protein[qf] for qf in quantfields}
    return qmap
