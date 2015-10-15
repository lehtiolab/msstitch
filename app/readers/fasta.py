from Bio import SeqIO


def get_proteins_for_db(fastafn):
    """Runs through fasta file and returns proteins accession nrs, sequences
    and evidence levels for storage in lookup DB. Duplicate accessions in
    fasta are accepted and removed by keeping only the last one.
    """
    objects = {}
    for record in parse_fasta(fastafn):
        objects[parse_protein_identifier(record)] = record
    return (((acc,) for acc in list(objects)),
            ((acc, str(record.seq)) for acc, record in objects.items()),
            ((acc, get_uniprot_evidence_level(record.description))
             for acc, record in objects.items()))


def generate_proteins_id(fastafn):
    for record in parse_fasta(fastafn):
        yield record.id


def get_proteins_descriptions(fastafn):
    for record in parse_fasta(fastafn):
        yield (record.id, record.description)


def get_proteins_genes(fastafn):
    """This returns a tuple of (protein, gene, HGNC symbol, description) from
    a passed file. If the file is FASTA from ENSEMBL or UniProt, only genes
    are given and symbol/description will be None. If the file is a ENSEMBL
    Biomart mapping file it tries to parse that and also return the rest"""
    with open(fastafn) as fp:
        firstline = next(fp).strip()
    if firstline[0] == '>':
        for record in parse_fasta(fastafn):
            rectype = get_record_type(record)
            yield (record.id, get_gene(record.description, rectype),
                   None, None)
    elif firstline.split('\t')[0] == 'Ensembl Gene ID':
        with open(fastafn) as fp:
            header = next(fp).strip().split('\t')
            ensg = header.index('Ensembl Gene ID')
            ensp = header.index('Ensembl Protein ID')
            symb = header.index('HGNC Symbol')
            desc = header.index('Description')
            for line in fp:
                line = line.strip().split('\t')
                yield (line[ensp], line[ensg], line[symb], line[desc])


def parse_protein_identifier(record):
    return record.id


def parse_fasta(fn):
    with open(fn) as fp:
        for record in SeqIO.parse(fp, 'fasta'):
            yield record


def get_record_type(record):
    if record.id.split('|')[0] == 'sp':
        return 'swiss'
    elif record.id[:4] == 'ENSP':
        return 'ensembl'
    else:
        raise RuntimeError('Cannot detect type of FASTA file. '
                           'Should be Uniprot or ENSEMBL')


def get_gene(description, rectype):
    splitter = {'ensembl': ':',
                'swiss': '='}[rectype]
    field = {'ensembl': 'gene',
             'swiss': 'GN'}[rectype]
    splitdesc = [x.split(splitter) for x in description.split()]
    return [x[1] for x in splitdesc if x[0] == field][0]


def has_evidence_levels(fastafn):
    fasta = parse_fasta(fastafn)
    record = next(fasta)
    if get_uniprot_evidence_level(record.description):
        return True
    return False


def get_uniprot_evidence_level(header):
    """Returns uniprot protein existence evidence level for a fasta header.
    Evidence levels are 1-5, but we return 5 - x since sorting still demands
    that higher is better."""
    header = header.split()
    for item in header:
        item = item.split('=')
        try:
            if item[0] == 'PE':
                return 5 - int(item[1])
        except IndexError:
            continue
    return -1
