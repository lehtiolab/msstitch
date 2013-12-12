from Bio import SeqIO
# import sqlite3


def get_searchspace(dbfns):
    """Given FASTA databases, proteins are trypsinized and resulting peptides
    stored in a database or dict for lookups"""
    allpepdict = {}
    for dbfn in dbfns:
        # allpeps = []
        protindex = SeqIO.index(dbfn, 'fasta')
        for acc in protindex:
            pepseqs = trypsinize(protindex[acc].seq)
            for pep in pepseqs:
                allpepdict[str(pep)] = 1
        ########## below here, we do if there is not enough memory on nodes.
        ## TEST FIRST WITH NORMAL DICT
            #allpeps.extend(pepseqs)
            #if len(allpeps)>10000: # more than x peps, then write to SQLite
                #write_peps_to_sqlite(allpeps)
                #allpeps = []
        # write remaining peps to sqlite
        #write_peps_to_sqlite(allpeps)
    return allpepdict

    
def write_peps_to_sqlite(peps):
    pass


def trypsinize(proseq):
    """Trypsinize a sequence. From Yafeng Zhu."""
    indice=[0]
    peptides=[]
    for i in range(0,len(proseq)-1):
        if proseq[i]=='K' and proseq[i+1]!='P':
            indice.append(i+1)
        elif proseq[i]=='R' and proseq[i+1]!='P':
            indice.append(i+1)

    indice.append(-1)                  
    for j in range(0,len(indice)-1):
        if indice[j]+1==indice[j+1]:
            peptides.append(proseq[indice[j]:indice[j+2]])
            peptides.append(proseq[indice[j-1]:indice[j+1]])
        else:
            peptides.append(proseq[indice[j]:indice[j+1]])
    return peptides


