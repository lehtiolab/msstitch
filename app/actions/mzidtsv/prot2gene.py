from app.readers import fasta as fastareader
from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata


def create_header(oldheader):
    p_ix = oldheader.index(mzidtsvdata.HEADER_PROTEIN) + 1
    return oldheader[:p_ix] + [mzidtsvdata.HEADER_GENE] + oldheader[p_ix:]


def add_genes_to_psm_table(psmfn, oldheader, fastafn):
    gpmap = get_protein_gene_map(fastafn)
    for psm in tsvreader.generate_tsv_psms(psmfn, oldheader):
        outpsm = {x: y for x, y in psm.items()}
        proteins = tsvreader.get_proteins_from_psm(psm)
        genes = get_genes(proteins, gpmap)
        genes = {gene: 1 for gene in genes.split(';')}
        outpsm[mzidtsvdata.HEADER_GENE] = ';'.join(genes)
        yield outpsm


def get_protein_gene_map(fastafn):
    gpmap = {}
    for protein, gene in fastareader.get_proteins_genes(fastafn):
        gpmap[protein] = gene
    return gpmap


def get_genes(proteins, gpmap):
    return ';'.join([gpmap[protein] for protein in proteins])
