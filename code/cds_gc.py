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

try:
    file = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <fasta file of coding regions>")



fa = rob.readFasta(file)
gc={}
total={}
for id in fa:
    pieces = id.split(" ")
    locus = pieces[1]
    if locus not in gc:
        gc[locus]=0
        total[locus]=0
    
    total[locus]+=len(fa[id])
    gc[locus]+=fa[id].lower().count('c')
    gc[locus]+=fa[id].lower().count('g')


print "Locus\tGC\tTotal\t%GC"

for locus in gc:
    print "\t".join([locus, str(gc[locus]), str(total[locus]), str(1.0*gc[locus]/total[locus])])



