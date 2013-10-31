"""
Converter script for Percolator output to create Fido-readable files.

By Jose Fernandez Navarro


usage (from galaxy):
percolator2fido.py $pout $graph
"""

import sys

from lxml import etree

class Hit(object):
    """Wrapper class"""
    def __init__(self, isdecoy=False,score=0.0,pep=0.0,q=0.0,p=0.0,name="",proteins=[]):
        self.isdecoy = isdecoy
        self.score = score
        self.pep = pep
        self.q = q
        self.p = p
        self.name = name
        self.proteins = proteins

class Perco2Fido(object):
    def __init__(self, args):
        infile = args[0]
        outfile = args[1]
        
        parser = etree.XMLParser(ns_clean=False, huge_tree=True)
        tree = etree.parse(infile, parser)
        elems = tree.getroot()  ##elems is what you pass to the readXXX functions
        
        peptides, pepunique = self.readPercolatorPeptides(elems)
        proteins = self.readPercolatorProteins(elems)
        self.writeFidoInput(peptides, proteins, outfile)
        
        
    def writeFidoInput(self, peptides,proteins,fidoOutputFile): #,fidoOutputFile2): # current released fido only outputs one file
    #   graph input looks like this:
    #    e EEEMPEPK
    #    r SW:TRP6_HUMAN
    #    r GP:AJ271067_1
    #    r GP:AJ271068_1
    #    p 0.9849
    #    e LLEIIQVR
    #    r SW:TRP6_HUMAN
    #    r GP:AJ271067_1
    #    r GP:AJ271068_1
    #    p 0.0
    #    e FAFNNKPNLEWNWK
    #    r gi|1574458|gb|AAC23247.1|
    #    p 0.9750

    #{ target1 , target2 , target3 , ... }
    #{ decoy1 , decoy2 , decoy3 , ... }

        f = open(fidoOutputFile, "w")
        for peptide in peptides:
            #if(not peptide.isdecoy):
            prots = peptide.proteins
            pepname = peptide.name
            pepprob = 1 - peptide.pep
            f.write("e " + pepname + "\n")
            for prot in prots:
                f.write("r " + prot + "\n")
            f.write("p " + str(pepprob) + "\n")    
        f.close()

######## implement this possibly with a newer fido version:
#        f = open(fidoOutputFile2, "w")
#        f.write("{ ")
#        for protein in [x for x in proteins if not x.isdecoy]:
#            f.write(printable(protein.name) + " , ")
#        f.write("}\n")
#        f.write("{ ")
#        for protein in [x for x in proteins if x.isdecoy]:
#            f.write(printable(protein.name) + " , ")
#        f.write("}\n")
#        f.close()



    def readPercolatorPSM(self, elems):
        ##process percolator Peptides,ptms and proteins
        percolatorPSMs = []
        percolatorPSMdict = dict()
        for elem in elems.iter("{http://per-colator.com/percolator_out/13}psm"):
            decoy = True   
            if (elem.get("{http://per-colator.com/percolator_out/13}decoy") == "false"):
                decoy = False
            score = elem.findtext("{http://per-colator.com/percolator_out/13}svm_score")
            pep = elem.findtext("{http://per-colator.com/percolator_out/13}pep")
            q = elem.findtext("{http://per-colator.com/percolator_out/13}q_value")
            p = elem.findtext("{http://per-colator.com/percolator_out/13}q_value")
            name = elem.findall("{http://per-colator.com/percolator_out/13}peptide_seq")[0].get("seq")
            percolatorPSMs.append(Hit(decoy, float(score), float(pep), float(q), float(p), str(name), []))
            if(not decoy):
                if(not percolatorPSMdict.has_key(name)):
                    percolatorPSMdict[name] = Hit(decoy, float(score), float(pep), float(q), float(p), str(name), [])
                elif(percolatorPSMdict[name].pep > float(pep)):
                    percolatorPSMdict[name].pep = float(pep)
                    percolatorPSMdict[name].score = float(score)
                    percolatorPSMdict[name].q = float(q)
                    percolatorPSMdict[name].p = float(p)
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
           
            return percolatorPSMs, percolatorPSMdict


    def readPercolatorPeptides(self, elems):
        percolatorPeptides = []
        peptideUnique = dict()
        for elem in elems.iter("{http://per-colator.com/percolator_out/13}peptide"):
            decoy = True   
            if (elem.get("{http://per-colator.com/percolator_out/13}decoy") == "false"):
                decoy = False
            score = elem.findtext("{http://per-colator.com/percolator_out/13}svm_score")
            pep = elem.findtext("{http://per-colator.com/percolator_out/13}pep")
            q = elem.findtext("{http://per-colator.com/percolator_out/13}q_value")
            p = elem.findtext("{http://per-colator.com/percolator_out/13}q_value")
            name = elem.get("{http://per-colator.com/percolator_out/13}peptide_id")
            proteins = elem.findall("{http://per-colator.com/percolator_out/13}protein_id")
            percolatorPeptides.append(Hit(decoy,float(score),float(pep),float(q),float(p),str(name),[x.text for x in proteins]))
            if(not decoy):
                if(not peptideUnique.has_key(str(name)) ):
                    peptideUnique[str(name)] = Hit(decoy,float(score),float(pep),float(q),float(p),str(name),[x.text for x in proteins])
                else:
                    print "Repeated peptide : " + name + " Found in Percolator files"
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
       
        return percolatorPeptides,peptideUnique


    def readPercolatorProteins(self, elems):
        percolatorProteins = []
        for elem in elems.iter("{http://per-colator.com/percolator_out/13}protein"):
            decoy = True   
            if (elem.get("{http://per-colator.com/percolator_out/13}decoy") == "false"):
                decoy = False
            pep = elem.findtext("{http://per-colator.com/percolator_out/13}pep")
            q = elem.findtext("{http://per-colator.com/percolator_out/13}q_value")
            p = elem.findtext("{http://per-colator.com/percolator_out/13}p_value")
            name = elem.get("{http://per-colator.com/percolator_out/13}protein_id")
            qmp = elem.findtext("{http://per-colator.com/percolator_out/13}q_value_emp")
            if (not qmp):
                qmp = 0
            percolatorProteins.append(Hit(decoy,float(qmp),float(pep),float(q),float(p),str(name),[]))
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
       
        return percolatorProteins

def main():
    convert = Perco2Fido(sys.argv[1:])
        
if __name__ == '__main__':
    main()

