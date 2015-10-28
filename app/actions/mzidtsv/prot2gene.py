from collections import OrderedDict

from app.readers import fasta as fastareader
from app.readers import tsv as tsvreader
from app.dataformats import mzidtsv as mzidtsvdata


def create_header(oldheader):
    p_ix = oldheader.index(mzidtsvdata.HEADER_PROTEIN) + 1
    newfields = [mzidtsvdata.HEADER_GENE, mzidtsvdata.HEADER_SYMBOL,
                 mzidtsvdata.HEADER_DESCRIPTION]
    return oldheader[:p_ix] + newfields + oldheader[p_ix:]


def add_genes_to_psm_table(psmfn, oldheader, pgdb):
    gpmap = pgdb.get_protein_gene_map()
    for psm in tsvreader.generate_tsv_psms(psmfn, oldheader):
        outpsm = {x: y for x, y in psm.items()}
        proteins = tsvreader.get_proteins_from_psm(psm)
        outpsm[mzidtsvdata.HEADER_GENE] = ';'.join(get_genes(proteins, gpmap))
        symbols = get_symbols(proteins, gpmap)
        desc = get_descriptions(proteins, gpmap)
        outpsm[mzidtsvdata.HEADER_SYMBOL] = ';'.join(symbols)
        outpsm[mzidtsvdata.HEADER_DESCRIPTION] = ';'.join(desc)
        yield outpsm


def get_mapped(proteins, gpmap, outtype):
    # FIXME multiple vals for a protein?
    outvals = [gpmap[protein][outtype] for protein in proteins]
    outvals = OrderedDict([(val, 1) for val in outvals]).keys()
    if None in outvals:
        return ['NA']
    else:
        return outvals


def get_genes(proteins, gpmap):
    return get_mapped(proteins, gpmap, 'gene')


def get_descriptions(proteins, gpmap):
    descriptions = get_mapped(proteins, gpmap, 'desc')
    return [x.replace('\t', ' ') for x in descriptions]


def get_symbols(proteins, gpmap):
    return get_mapped(proteins, gpmap, 'symbol')
