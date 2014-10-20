from Bio import SeqIO


def get_proteins_for_db(fastafn, evidence_levels=False):
    """Convenience function to get 3 iterators of fasta file"""
    if evidence_levels:
        evidence_levels = get_evidence_iter(fastafn)
    return (get_protein_acc_iter(fastafn), get_sequence_iter(fastafn),
            evidence_levels)


def get_protein_acc_iter(fastafn):
    """Returns iterator with protein accessions"""
    for record in parse_fasta(fastafn):
        yield record.id


def get_sequence_iter(fastafn):
    """Returns iterator with protein accessions, sequences in tuple"""
    for record in parse_fasta(fastafn):
        yield (record.id, record.seq)


def get_evidence_iter(fastafn):
    """Returns iterator with protein accessions, evidence levels in tuple"""
    for record in parse_fasta(fastafn):
        yield (record.id, get_uniprot_evidence_level(record.description))


def parse_fasta(fn):
    with open(fn) as fp:
        return SeqIO.parse(fp, 'fasta')


def has_evidence_levels(fastafn):
    fasta = parse_fasta(fastafn)
    record = next(fasta)
    if get_uniprot_evidence_level(record.description):
        return True
    return False


def get_uniprot_evidence_level(header):
    header = header.split()
    for item in header:
        item = item.split('=')
        try:
            if item[0] == 'PE':
                return item[1]
        except IndexError:
            continue
    return False
