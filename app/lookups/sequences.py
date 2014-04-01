from Bio import SeqIO
from app import sqlite


def create_searchspace(dbfns):
    """Given FASTA databases, proteins are trypsinized and resulting peptides
    stored in a database or dict for lookups"""
    lookup = sqlite.SearchSpaceDB()
    lookup.create_searchspacedb()

    for dbfn in dbfns:
        allpeps = []
        protindex = SeqIO.index(dbfn, 'fasta')
        for acc in protindex:
            pepseqs = trypsinize(protindex[acc].seq)
            # Exchange all leucines to isoleucines because MS can't differ
            pepseqs = [(str(pep).replace('L', 'I'),) for pep in pepseqs]
            allpeps.extend(pepseqs)
            if len(allpeps) > 1000000:  # more than x peps, write to SQLite
                lookup.write_peps(allpeps)
                allpeps = []
        # write remaining peps to sqlite
        lookup.write_peps(allpeps)
    lookup.index_peps()
    lookup.close_connection()
    return lookup.fn


def trypsinize(proseq):
    """Trypsinize a sequence. From Yafeng Zhu. Returns fully
    tryptic, and overlapping peptides if there are a sequence of
    tryptic residues"""
    indices = [0]
    peptides = []
    for aa, i in enumerate(range(0, len(proseq) - 1)):
        nextaa = proseq[i + 1]
        if aa in ['K', 'P'] and nextaa != 'P':
            indices.append(i + 1)

    indices.append(-1)
    for j in range(0, len(indices) - 1):
        if indices[j] + 1 == indices[j + 1]:
            peptides.append(proseq[indices[j]:indices[j + 2]])
            peptides.append(proseq[indices[j - 1]:indices[j + 1]])
        elif indices[j + 1] == -1:
            peptides.append(proseq[indices[j]:])
        else:
            peptides.append(proseq[indices[j]:indices[j + 1]])
    return peptides
