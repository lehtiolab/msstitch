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


def create_searchspace(lookup, fastafn, proline_cut=False, reverse_seqs=True,
        do_trypsinize=True, fully_tryptic=False):
    """Given a FASTA database, proteins are trypsinized and resulting peptides
    stored in a database or dict for lookups"""
    allpeps = []
    for record in SeqIO.parse(fastafn, 'fasta'):
        if do_trypsinize:
            pepseqs = trypsinize(record.seq, proline_cut, fully_tryptic=fully_tryptic)
        else:
            pepseqs = [record.seq]
        # Exchange all leucines to isoleucines because MS can't differ
        pepseqs = [(str(pep).replace('L', 'I'),) for pep in pepseqs]
        allpeps.extend(pepseqs)
        if len(allpeps) > 1000000:  # more than x peps, write to SQLite
            lookup.write_peps(allpeps, reverse_seqs)
            allpeps = []
    # write remaining peps to sqlite
    lookup.write_peps(allpeps, reverse_seqs)
    lookup.index_peps(reverse_seqs)
    #lookup.close_connection()


def trypsinize(proseq, proline_cut=False, fully_tryptic=False):
    # TODO add cysteine to non cut options
    """Trypsinize a protein sequence. Returns a list of peptides.
    Peptides include both cut and non-cut when P is behind a tryptic
    residue. Multiple consequent tryptic residues are treated as follows:
    PEPKKKTIDE - [PEPK, PEPKK, PEPKKK, KKTIDE, KTIDE, TIDE, K, K, KK ]
    When fully_tryptic = True, the peptide would yield: [PEPK, K, K, TIDE]
    """
    outpeps = []
    currentpeps = ['']
    trypres = set(['K', 'R'])
    noncutters = set()
    if not proline_cut:
        noncutters.add('P')
    for i, aa in enumerate(proseq):
        currentpeps = ['{0}{1}'.format(x, aa) for x in currentpeps]
        if i == len(proseq) - 1:
            continue
        if aa in trypres and proseq[i + 1] not in noncutters:
            outpeps.extend(currentpeps)  # do actual cut by storing peptides
            if fully_tryptic:
                currentpeps = ['']
            elif proseq[i + 1] in trypres.union('P'):
                # add new peptide to list if we are also to run on
                currentpeps.append('')
            elif trypres.issuperset(currentpeps[-1]):
                currentpeps = [x for x in currentpeps if trypres.issuperset(x)]
                currentpeps.append('')
            else:
                currentpeps = ['']

    if currentpeps != ['']:
        outpeps.extend(currentpeps)
    return outpeps


def tryp_rev(seq, lookup, do_trypsinize):
    if do_trypsinize:
        segments = trypsinize(seq, fully_tryptic=True)
    else:
        segments = [str(seq.seq)]
    final_seq = []
    for s in segments :
        if len(s) > 1 :
            if s[-1] in ['R', 'K']:
                new_s = '{}{}'.format(s[:-1][::-1], s[-1])
            else:
                new_s = s[::-1]
            shufflecount = 0
            while lookup and lookup.check_seq_exists(new_s.replace('L', 'I'), amount_ntermwildcards=0) and shufflecount < 100:
                nterm = list(new_s[:-1])
                shuffle(nterm)
                new_s = '{}{}'.format(''.join(nterm), new_s[-1])
                shufflecount += 1
            if shufflecount < 10: # max 10 shuffles, else discard decoy
                final_seq.append(new_s)
        else :
            new_s = s
            final_seq.append(s)
    if final_seq:
        seq.seq = Seq(''.join(final_seq))
        seq.id = 'decoy_{}'.format(seq.name)
    else:
        seq = False
    return seq


def prot_rev(seq, lookup):
    seq.id = 'decoy_{}'.format(seq.name)
    seq = seq[::-1]
    return seq


def create_decoy_fa(fastafn, method, lookup, is_trypsinized):
    outfasta = SeqIO.parse(fastafn, 'fasta')
    if method == 'prot_rev':
        outfasta = (prot_rev(x, lookup) for x in outfasta)
    if method == 'tryp_rev':
        outfasta = (tryp_rev(x, lookup, is_trypsinized) for x in outfasta)
    return (x for x in outfasta if x) # do not yield empty records
