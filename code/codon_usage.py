''' 
Generate a set of vectors of codon usage for different organisms. We can then compare them by Euclidean distance.

We accept a fasta file typically created by ~/bioinformatics/ncbi/combine_gbff_fna.py which has the following format of the fasta line
>LOCUS_TAG LOCUS LOCATION

The LOCUS is the genome identifier we need to use

'''

import os
import sys
import numpy as np
sys.path.append('/home3/redwards/bioinformatics/Modules')
import rob
import robseq

# this is a hash of all possible codons
codons = robseq.geneticCode().keys()
codons.sort()

try:
    file = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <fasta file of coding regions>")



fa = rob.readFasta(file)
count={}
cds={}
for id in fa:
    if ((1.0 * len(fa[id])) / 3) != len(fa[id])/3:
        sys.stderr.write("Sequence " + id + " does not appear to be a multiple of 3 nucleotides. Skipped\n")
        continue
    pieces = id.split(" ")
    locus = pieces[1]
    if locus not in count:
        count[locus]=0
        cds[locus]={}
        for codon in codons:
            cds[locus][codon]=0
    p=0
    while p<len(fa[id]):
        c = fa[id][p:p+3].upper()
        p += 3
        if c not in codons:
            # this is most likely a sequence ambiguity that does not make a real amino acid
            # ignore for now 
            continue
        cds[locus][c] += 1
        count[locus] += 1

print "Locus",
for codon in codons:
    print "\t" + codon,
print 

for locus in count:
    print locus,
    for codon in codons:
        print "\t" + str(1.0 * cds[locus][codon] / count[locus]),
    print



