from Bio import SeqIO


def get_proteins_for_db(fastafn, evidence_levels=False):
    """Runs through fasta file and returns proteins accession nrs, sequences
    and evidence levels for storage in lookup DB. Duplicate accessions in
    fasta are accepted and removed by keeping only the last one.
    """
    objects = {}
    for record in parse_fasta(fastafn):
        objects[parse_protein_identifier(record)] = record
    if evidence_levels:
        return (((acc,) for acc in list(objects)),
                ((acc, str(record.seq)) for acc, record in objects.items()),
                ((acc, get_uniprot_evidence_level(record.description))
                 for acc, record in objects.items()))
    else:
        return (((acc,) for acc in list(objects)),
                ((acc, str(record.seq)) for acc, record in objects.items()),
                False)


def get_proteins_descriptions(fastafn):
    for record in parse_fasta(fastafn):
        yield (record.id, record.desc)


def parse_protein_identifier(record):
    return record.id


def parse_fasta(fn):
    with open(fn) as fp:
        for record in SeqIO.parse(fp, 'fasta'):
            yield record


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
                return 5 - item[1]
        except IndexError:
            continue
    return False
