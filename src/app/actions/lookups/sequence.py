from Bio import SeqIO
from collections import defaultdict

from app.actions.sequence import trypsinize

PROTEIN_STORE_CHUNK_SIZE = 100000


def create_searchspace_wholeproteins(lookup, fastafn, minpeplen):
    fasta = SeqIO.parse(fastafn, 'fasta')
    storeseqs = {}
    peptotal, pepsubtotal = 0, 0
    for record in fasta:
        prot_id, protseq = record.id, str(record.seq)
        storeseqs[prot_id] = {'seq': protseq, 'peps': []}
        isoseq = protseq.replace('L', 'I') 
        for pos in range(0, len(protseq) - minpeplen + 1):
            possible_pep = isoseq[pos:pos + minpeplen]
            pepsubtotal += 1
            storeseqs[prot_id]['peps'].append((possible_pep, pos))
        if pepsubtotal > PROTEIN_STORE_CHUNK_SIZE:
            lookup.store_pep_proteins(storeseqs)
            peptotal += pepsubtotal
            pepsubtotal = 0
            storeseqs = {}
    peptotal += pepsubtotal
    lookup.store_pep_proteins(storeseqs)
    lookup.index_proteins()
    print(f'Stored {peptotal} peptides')


def process_sequence(seq, do_trypsinize, proline_cut, miss_cleavage, minlen, reverse, ntermmloss): 
    if do_trypsinize:
        pepseqs = trypsinize(seq, proline_cut, miss_cleavage=miss_cleavage, split_stop_codons=True,
                opt_nt_meth_loss=ntermmloss)
    else:
        pepseqs = [(seq)]
    # Exchange all leucines to isoleucines because MS can't differ
    pepseqs = (pep.replace('L', 'I') for pep in pepseqs)
    if minlen:
        pepseqs = (pep for pep in pepseqs if len(pep) >= minlen)
    if reverse:
        pepseqs = (pep[::-1] for pep in pepseqs)
    return pepseqs


def create_searchspace(lookup, infile, minlen, proline_cut=False, reverse_seqs=True,
        do_trypsinize=True, miss_cleavage=False, ntermmloss=False):
    """Given a FASTA database, proteins are trypsinized and resulting peptides
    stored in a database or dict for lookups"""
    with open(infile) as fp:
        ftype = 'fasta' if fp.read(1) == '>' else 'txt'
        fp.seek(0)
        if ftype == 'fasta':
            records = (str(x.seq) for x in SeqIO.parse(fp, 'fasta'))
        elif ftype == 'txt':
            records = (x.strip('\n') for x in fp)
        lookup.write_peps((pep,) for x in records for pep in process_sequence(x, do_trypsinize,
            proline_cut, miss_cleavage, minlen, reverse_seqs, ntermmloss))
    lookup.index_peps(reverse_seqs)


def create_searchspace_map_accessions(lookup, infile, minlen, proline_cut=False, reverse_seqs=True,
        do_trypsinize=True, miss_cleavage=False, ntermmloss=False):
    """Given a FASTA database, proteins are trypsinized and resulting peptides
    stored in a database or dict for lookups, but this time with mapping protein accessions
    also stored, so we can find them later"""
    allpeps = defaultdict(list)
    peptotal, pepsubtotal = 0, 0
    with open(infile) as fp:
        for rec in SeqIO.parse(fp, 'fasta'):
            peps = process_sequence(rec.seq, do_trypsinize, proline_cut, miss_cleavage, minlen,
                reverse_seqs, ntermmloss)
            allpeps[rec.id].extend(peps)
    lookup.store_tryp_peps_mapped(allpeps) 
    lookup.index_peps_mapped(reverse_seqs)
