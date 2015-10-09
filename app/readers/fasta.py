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
    for record in parse_fasta(fastafn):
        rectype = get_record_type(record)
        yield (record.id, get_gene(record.description, rectype))


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
        raise RuntimeError('Cannot detect type of FASTA file. Should be Uniprot or ENSEMBL')


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
