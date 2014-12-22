'''
Choose a set of genomes with highly expressed proteins, and write just those ORFs to a file so we can build a codon adaptation index from them.

The files you probably want are:
    bacteria_taxid.txt 
    highly_expressed_proteins.completegenomes.blastp
    /lustre/usr/data/NCBI/RefSeq/bacteria/refseq_orfs.faa
'''

import sys
import os
import re

try:
    ncF = sys.argv[1]
    bhF = sys.argv[2]
    rfH = sys.argv[3]
    oF  = sys.argv[4]
    nOrfs = int(sys.argv[5])
except:
    sys.exit(sys.argv[0] + " <List of genomes we want to use> <blast hits file> <refseq ORFs file> <ORF output file> <minimum number of ORFs that must be present to include a genome>\nGenomes should be NC numbers, no versions, one per line (you can use the taxonomy file bacteria_taxid.txt.\nBlast hits file should be blastp of protein sequences against highly expressed sequences\nRefSeq ORFs file should have the same ids as the proteins\n")

# find just the set of genomes we are looking for
want={}
with open(ncF, 'r') as ncin:
    for l in ncin:
        p=l.strip().split('\t')
        want[p[0]]=1

# identify the highly expressed proteins by similarity to a known set
he={}
with open(bhF, 'r') as bhin:
    for l in bhin:
        p=l.strip().split('\t')
        he[p[1]]=1

# finally, find the ORFs we are looking for
count={}
with open(rfH, 'r') as rin:
    for l in rin:
        if l.startswith('>'):
            m = re.findall('>([\w\.\-]+)\s+\[(\w+)\]\s+\[(.*)\]', l)
            if m == []:
                sys.stderr.write("ERROR parsing " + l)
                continue
            gene, genome, locus = m[0]
            if gene in he and genome in want:
                count[genome]=count.get(genome, 0) + 1

p=0
with open(rfH, 'r') as rin:
    with open(oF, 'w') as out:
        for l in rin:
            if l.startswith('>'):
                p=0
                m = re.findall('>([\w\.\-]+)\s+\[(\w+)\]\s+\[(.*)\]', l)
                if m == []:
                    continue
                gene, genome, locus = m[0]
                if gene in he and genome in want and count[genome] >= nOrfs:
                    p=1
            if p > 0:
                out.write(l)


for g in count:
    if count[g] > nOrfs:
        print g + "\t" + str(count[g])


