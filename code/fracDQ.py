'''

Implement the Frac_Q and Frac_D methods from http://cge.cbs.dtu.dk/services/HostPhinder/.

We start with a directory of text files, where each genome has a list of kmers that are found in that genome. Then we need to know the number shared. That phage then predicts the hos.

'''

import sys,os
from phage import Phage
phage=Phage()
phages=phage.phageTaxonomy()

try:
    kmerD=sys.argv[1]
    outf=sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <15-mer host directory> <output file>")

# read in all the 15mers
kmer={}
allk={}
for p in phages:
    kmer[p]={}
    if not os.path.exists(os.path.join(kmerD, p + ".tsv")):
        sys.stderr.write("No kmer file for " + p + "\n")

    with open(os.path.join(kmerD, p + ".tsv"), 'r') as fin:
        for l in fin:
            part=l.split("\t")
            kmer[p][part[0]]=1
            if part[0] not in allk:
                allk[part[0]] = {}
            allk[part[0]][p]=1

phgs=phages.keys()
phgs.sort()
with open(outf, 'w') as out:
    out.write("#Phage 1\tPhage 2\tFrac_Q\tFrac_D\n")
    for i in xrange(len(phgs)-1):
        for j in xrange(i+1, len(phgs)):
            # how many 15mers do these share?
            share=0
            for k in kmer[phgs[i]]:
                if phgs[j] in allk[k]:
                    share+=1
            # frac_Q assumes phage 1 is query and phage 2 is database
            out.write("\t".join([phgs[i], phgs[j]]))
            out.write("\t" + str(1.0*share/len(kmer[phgs[i]])))
            out.write("\t" + str(1.0*share/len(kmer[phgs[j]])))
            out.write("\n")
            # now print out for the other way around
            out.write("\t".join([phgs[j], phgs[i]]))
            out.write("\t" + str(1.0*share/len(kmer[phgs[j]])))
            out.write("\t" + str(1.0*share/len(kmer[phgs[i]])))
            out.write("\n")






