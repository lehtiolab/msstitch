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
    """Trypsinize a sequence. From Yafeng Zhu."""
    indice = [0]
    peptides = []
    for i in range(0, len(proseq) - 1):
        if proseq[i] == 'K' and proseq[i + 1] != 'P':
            indice.append(i + 1)
        elif proseq[i] == 'R' and proseq[i + 1] != 'P':
            indice.append(i + 1)

    indice.append(-1)
    for j in range(0, len(indice) - 1):
        if indice[j] + 1 == indice[j + 1]:
            peptides.append(proseq[indice[j]:indice[j + 2]])
            peptides.append(proseq[indice[j - 1]:indice[j + 1]])
        else:
            peptides.append(proseq[indice[j]:indice[j + 1]])
    return peptides
