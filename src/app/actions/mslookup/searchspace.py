from Bio import SeqIO
from Bio.Seq import Seq
from random import shuffle
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


def create_searchspace(lookup, fastafn, minlen, proline_cut=False, reverse_seqs=True,
        do_trypsinize=True, miss_cleavage=False):
    """Given a FASTA database, proteins are trypsinized and resulting peptides
    stored in a database or dict for lookups"""
    allpeps = []
    for record in SeqIO.parse(fastafn, 'fasta'):
        if do_trypsinize:
            pepseqs = trypsinize(record.seq, proline_cut, miss_cleavage=miss_cleavage)
        else:
            pepseqs = [record.seq]
        # Exchange all leucines to isoleucines because MS can't differ
        pepseqs = [(str(pep).replace('L', 'I'),) for pep in pepseqs]
        if minlen:
            pepseqs = [pep for pep in pepseqs if len(pep[0]) >= minlen]
        allpeps.extend(pepseqs)
        if len(allpeps) > 1000000:  # more than x peps, write to SQLite
            lookup.write_peps(allpeps, reverse_seqs)
            allpeps = []
    # write remaining peps to sqlite
    lookup.write_peps(allpeps, reverse_seqs)
    lookup.index_peps(reverse_seqs)
    #lookup.close_connection()


def trypsinize(proseq, proline_cut=False, miss_cleavage=0):
    # TODO add cysteine to non cut options
    """Trypsinize a protein sequence. Returns a list of peptides.
    Peptides include both cut and non-cut when P is behind a tryptic
    residue. Number of missed cleavages can be specified.
    """
    #if not set(str(proseq.seq[1:-1])).intersection('RK'):
    #    return [str(proseq.seq)]
    outpeps = []
    trypres = set(['K', 'R'])
    noncutters = set()
    if not proline_cut:
        noncutters.add('P')
    peptide = ''
    for i, aa in enumerate(proseq):
        peptide += aa
        if i == len(proseq) - 1:
            outpeps.append(peptide)
        elif aa in trypres and proseq[i + 1] not in noncutters:
            outpeps.append(peptide)  # do actual cut by storing peptides
            peptide = ''
        
    ft_pep_amount = len(outpeps)
    for i in range(1, miss_cleavage + 1):
        for j, pep in enumerate([x for x in outpeps]):
            if j + i < ft_pep_amount:
                outpeps.append(''.join([x for x in (outpeps[j:j + i + 1])]))
    return outpeps


def tryp_rev(seq, lookup, do_trypsinize, miss_cleavage, minlen):
    if do_trypsinize:
        segments = trypsinize(seq, miss_cleavage=miss_cleavage)
    else:
        segments = [str(seq.seq)]
    final_seq = {}
    decoy_segs = {}
    for i, s in enumerate(segments):
        if len(s) > 1 :
            if s[-1] in ['R', 'K']:
                decoy_segs[i] = '{}{}'.format(s[:-1][::-1], s[-1])
            else:
                decoy_segs[i] = s[::-1]
        else:
            decoy_segs[i] = s
    if lookup is not None:
        shufflecount = 0
        targets, tests = True, {k:v for k,v in decoy_segs.items()}
        if minlen:
            tests = {k:v for k,v in tests.items() if len(v) >= minlen}
        while targets and shufflecount < 10:
            targets = lookup.get_multi_seq(list(tests.values()))
            shuffled = 0
            for i,s in [(k,v) for k,v in tests.items()]:
                if s not in targets:
                    decoy_segs[i] = tests.pop(i)
                else:
                    nterm = list(s[:-1])
                    shuffle(nterm)
                    tests[i] = '{}{}'.format(''.join(nterm), s[-1])
                    shuffled = 1
            shufflecount += shuffled
        # discard max shuffled peptides
        if shufflecount >= 10:
            decoy_segs.update({i: '' for i in tests.keys()})
    if decoy_segs:
        seq.seq = Seq(''.join([decoy_segs[i] for i in range(0, len(decoy_segs))]))
        seq.id = 'decoy_{}'.format(seq.name)
    else:
        seq = False
    return seq


def prot_rev(seq, lookup):
    seq.id = 'decoy_{}'.format(seq.name)
    seq = seq[::-1]
    return seq


def create_decoy_fa(fastafn, method, lookup, is_trypsinized, miss_cleavage, minlen):
    outfasta = SeqIO.parse(fastafn, 'fasta')
    if method == 'prot_rev':
        outfasta = (prot_rev(x, lookup) for x in outfasta)
    if method == 'tryp_rev':
        outfasta = (tryp_rev(x, lookup, is_trypsinized, miss_cleavage, minlen) for x in outfasta)
    return (x for x in outfasta if x) # do not yield empty records
