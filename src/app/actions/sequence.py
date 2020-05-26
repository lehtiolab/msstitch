from Bio import SeqIO
from Bio.Seq import Seq
from random import shuffle


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
