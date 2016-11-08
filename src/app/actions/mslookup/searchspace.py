from Bio import SeqIO
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


def create_searchspace(lookup, fastafn, proline_cut=False,
                       reverse_seqs=True, do_trypsinize=True):
    """Given a FASTA database, proteins are trypsinized and resulting peptides
    stored in a database or dict for lookups"""
    allpeps = []
    for record in SeqIO.parse(fastafn, 'fasta'):
        if do_trypsinize:
            pepseqs = trypsinize(record.seq, proline_cut)
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
    lookup.close_connection()


def trypsinize(proseq, proline_cut=False):
    # TODO add cysteine to non cut options, use enums
    """Trypsinize a protein sequence. Returns a list of peptides.
    Peptides include both cut and non-cut when P is behind a tryptic
    residue. Multiple consequent tryptic residues are treated as follows:
    PEPKKKTIDE - [PEPK, PEPKK, PEPKKK, KKTIDE, KTIDE, TIDE, K, K, KK ]
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
            if proseq[i + 1] in trypres.union('P'):
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
