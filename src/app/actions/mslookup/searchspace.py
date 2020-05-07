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
    #lookup.close_connection()


def create_trypsinized(proteins, proline_cut=False, miss_cleavage=0, minlen=0):
    for proseq in proteins:
        trypnr = 0
        outpeps = trypsinize(str(proseq.seq), proline_cut, miss_cleavage)
        if minlen:
            outpeps = [pep for pep in outpeps if len(pep) >= minlen]
        for pep in outpeps:
            trypnr += 1
            trypseq = SeqIO.SeqRecord(Seq(pep), id='{}_{}'.format(proseq.id, trypnr), name=proseq.name, description='')
            yield trypseq


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


def tryp_rev(seq, lookup, do_trypsinize, miss_cleavage, minlen, max_shuffle):
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
    nr_decoymatching = 0
    if minlen:
        decoy_segs = {k: (v if len(v) >=minlen else '') for k,v in decoy_segs.items() }
    nr_peptides = len([x for x in decoy_segs.values() if x != ''])
    if lookup is not None:
        shufflecount = 0
        targets, tests = True, {k: (v, 0, v) for k,v in decoy_segs.items()}
        while targets:
            targets = lookup.get_multi_seq([x[0] for x in tests.values()])
            for i, (s, shufcount, origdecoy) in [(k,v) for k,v in tests.items()]: # list comprehension to not have dict change during iteration
                if s not in targets:
                    decoy_segs[i] = tests.pop(i)[0]
                elif shufcount < max_shuffle:
                    nterm = list(s[:-1])
                    shuffle(nterm)
                    tests[i] = ('{}{}'.format(''.join(nterm), s[-1]), shufcount + 1, origdecoy)
                else:
                    decoy_segs[i] = origdecoy
                    tests.pop(i)
                    nr_decoymatching += 1
    if set(decoy_segs.values()) != {''}:
        seq.seq = Seq(''.join([decoy_segs[i] for i in range(0, len(decoy_segs))]))
        seq.id = 'decoy_{}'.format(seq.name)
        seq.description = 'decoy_{}'.format(seq.description)
    else:
        seq = False
    return seq, nr_decoymatching, nr_peptides


def prot_rev(seq):
    seq.id = 'decoy_{}'.format(seq.name)
    seq = seq[::-1]
    return seq


def create_decoy_fa(fastafn, method, lookup, is_trypsinized, miss_cleavage, minlen, max_shuffle):
    outfasta = SeqIO.parse(fastafn, 'fasta')
    if method == 'prot_rev':
        for seq in outfasta:
            yield prot_rev(seq)
    if method == 'tryp_rev':
        decoymatching, nr_peptides = 0, 0
        outtryp = (tryp_rev(x, lookup, is_trypsinized, miss_cleavage, minlen, 
            max_shuffle) for x in outfasta)
        for seq, nr_decoymatching, nr_pep in outtryp:
            decoymatching += nr_decoymatching
            nr_peptides += nr_pep
            if seq:
                yield seq
        if decoymatching:
            print('Unable to shuffle {} decoys (of a total of {} decoy peptides) '
                    'that matched target DB (retained non-shuffled)'.format(decoymatching, nr_peptides))
