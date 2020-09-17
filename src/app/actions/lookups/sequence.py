from Bio import SeqIO

from app.actions.sequence import trypsinize

PROTEIN_STORE_CHUNK_SIZE = 100000


def create_searchspace_wholeproteins(lookup, fastafn, minpeplen):
    fasta = SeqIO.parse(fastafn, 'fasta')
    prots = {str(prot.seq).replace('L', 'I'): prot.id for prot in fasta}
    storeseqs = []
    peptotal = 0
    for protseq, prot_id in prots.items():
        for pos in range(0, len(protseq) - minpeplen + 1):
            possible_pep = protseq[pos:pos + minpeplen]
            peptotal += 1
            storeseqs.append((possible_pep, prot_id, pos))
        if len(storeseqs) > PROTEIN_STORE_CHUNK_SIZE:
            lookup.store_pep_proteins(storeseqs)
            storeseqs = []
    lookup.store_pep_proteins(storeseqs)
    lookup.index_proteins()
    print('Stored {} peptides from {} proteins (reduced FASTA to remove '
          'duplicate sequences)'.format(peptotal, len(prots)))


def create_searchspace(lookup, infile, minlen, proline_cut=False, reverse_seqs=True,
        do_trypsinize=True, miss_cleavage=False):
    """Given a FASTA database, proteins are trypsinized and resulting peptides
    stored in a database or dict for lookups"""
    allpeps = []
    def treat_pep(seq, do_trypsinize, proline_cut, miss_cleavage, minlen): 
        if do_trypsinize:
            pepseqs = trypsinize(seq, proline_cut, miss_cleavage=miss_cleavage)
        else:
            pepseqs = [seq]
        # Exchange all leucines to isoleucines because MS can't differ
        pepseqs = [(pep.replace('L', 'I'),) for pep in pepseqs]
        if minlen:
            pepseqs = [pep for pep in pepseqs if len(pep[0]) >= minlen]
        return pepseqs
    with open(infile) as fp:
        ftype = 'fasta' if fp.read(1) == '>' else 'txt'
        fp.seek(0)
        if ftype == 'fasta':
            records = (str(x.seq) for x in SeqIO.parse(fp, 'fasta'))
        elif ftype == 'txt':
            records = (x.strip('\n') for x in fp)
        lookup.write_peps((pep for x in records for pep in treat_pep(x, do_trypsinize, proline_cut, miss_cleavage, minlen)), reverse_seqs)
    lookup.index_peps(reverse_seqs)
