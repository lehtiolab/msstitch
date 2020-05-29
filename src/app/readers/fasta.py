from Bio import SeqIO


def get_proteins_for_db(fastafn, fastadelim, genefield):
    """Runs through fasta file and returns proteins accession nrs, sequences
    and evidence levels for storage in lookup DB. Duplicate accessions in
    fasta are accepted and removed by keeping only the last one.
    """
    records = {acc: (rec, get_record_type(rec)) for acc, rec in 
            SeqIO.index(fastafn, 'fasta').items()}
    proteins = ((x,) for x in records.keys())
    sequences = ((acc, str(rec.seq)) for acc, (rec, rtype) in records.items())
    desc = ((acc, get_description(rec, rtype)) for acc, (rec, rtype) in records.items() if rtype)
    evid = ((acc, get_uniprot_evidence_level(rec, rtype)) for acc, (rec, rtype) in 
        records.items())
    ensgs = [(get_ensg(rec), acc) for acc, (rec, rtype) in records.items()
            if rtype == 'ensembl']
    def sym_out():
        symbols = ((get_symbol(rec, rtype, fastadelim, genefield), acc) for 
                acc, (rec, rtype) in records.items() if rtype)
        othergene = ((get_other_gene(rec, fastadelim, genefield), acc) for acc, (rec, rtype) in records.items()
                if not rtype and fastadelim and fastadelim in rec.description)
        yield from symbols
        yield from othergene
    return proteins, sequences, desc, evid, ensgs, [x for x in sym_out()]


def parse_fasta(fn):
    with open(fn) as fp:
        for record in SeqIO.parse(fp, 'fasta'):
            yield record


def get_record_type(record):
    dmod = get_decoy_mod_string(record.id)
    test_name = record.id
    if dmod is not None:
        test_name = record.id.replace(dmod, '')
    if test_name.split('|')[0] in ['sp', 'tr']:
        return 'swiss'
    elif test_name[:3] == 'ENS':
        return 'ensembl'
    else:
        return False


def get_decoy_mod_string(protein):
    mods = ['tryp_reverse', 'reverse', 'decoy', 'random', 'shuffle']
    for mod in mods:
        if mod in protein:
            if protein.endswith('_{}'.format(mod)):
                return '_{}'.format(mod)
            elif protein.endswith('{}'.format(mod)):
                return mod
            elif protein.startswith('{}_'.format(mod)):
                return '{}_'.format(mod)
            elif protein.startswith('{}'.format(mod)):
                return mod


def get_description(record, rectype):
    if rectype == 'ensembl':
        desc_spl = [x.split(':') for x in record.description.split()]
        try:
            descix = [ix for ix, x in enumerate(desc_spl) if x[0] == 'description'][0]
        except IndexError:
            return 'NA'
        desc = ' '.join([':'.join(x) for x in desc_spl[descix:]])[12:]
        return desc
    elif rectype == 'swiss':
        desc = []
        for part in record.description.split()[1:]:
            if len(part.split('=')) > 1:
                break
            desc.append(part)
        return ' '.join(desc)


def get_other_gene(record, fastadelim, genefield):
    return record.description.split(fastadelim)[genefield]


def get_genes_pickfdr(fastafn, outputtype, fastadelim, genefield):
    """Called by protein FDR module for both ENSG and e.g. Uniprot"""
    for rec in parse_fasta(fastafn):
        rtype = get_record_type(rec)
        if rtype == 'ensembl' and outputtype == 'ensg':
            yield get_ensg(rec)
        elif outputtype == 'genename':
            yield get_symbol(rec, rtype, fastadelim, genefield)


def get_ensg(record):
    fields = [x.split(':') for x in record.description.split()]
    try:
        return [x[1] for x in fields if x[0] == 'gene' and len(x) == 2][0]
    except IndexError:
        raise RuntimeError('ENSEMBL detected but cannot find gene ENSG in fasta')


def get_symbol(record, rectype, fastadelim, genefield):
    if rectype == 'ensembl':
        fields = [x.split(':') for x in record.description.split()]
        sym = [x[1] for x in fields if x[0] == 'gene_symbol' and len(x) == 2]
    elif rectype == 'swiss':
        fields = [x.split('=') for x in record.description.split()]
        sym = [x[1] for x in fields if x[0] == 'GN' and len(x) == 2]
    elif fastadelim and fastadelim in record.description and genefield:
        return record.description.split(fastadelim)[genefield]
    else:
        return 'NA'
    try:
        return sym[0]
    except IndexError:
        return 'NA'


def get_uniprot_evidence_level(record, rtype):
    """Returns uniprot protein existence evidence level for a fasta header.
    Evidence levels are 1-5, but we return 5 - x since sorting still demands
    that higher is better."""
    if rtype != 'swiss':
        return -1
    for item in record.description.split():
        item = item.split('=')
        try:
            if item[0] == 'PE' and len(item) == 2:
                return 5 - int(item[1])
        except IndexError:
            continue
    return -1
