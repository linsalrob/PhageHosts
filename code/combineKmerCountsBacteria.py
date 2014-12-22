#!/usr/bin/python

'''Read the kmer counts from a series of files and plot a PCA'''

## start with reading the files and adding the kmer counts to the organism names
import sys
sys.path.append('/home3/redwards/bioinformatics/phage_host')
from phage import Phage
import re
import os

phage = Phage()
## we only choose those hosts with 5 phages that infect them.
#host = phage.phageHost()
bacteria = phage.completeBacteria()


try:
    dir  = sys.argv[1]
    outf = sys.argv[2]
    merWanted = sys.argv[3]
except:
    sys.exit(sys.argv[0] + " <directory of kmer counts> <file to write output table to> <kmer size>")

count={}
allkmers={}
organismId={}
for file in os.listdir(dir):
    # our file names look kmers/165.kmers/NC_008025.1.3kmer.tsv
    match = re.match('(.*)\.\d+\.(\d+)kmer.tsv', file)
    if match ==None:
        sys.stderr.write("Did not find an id in " + file + " Skipped\n")
        continue

    id = match.group(1)
    mer = match.group(2)
    if mer != merWanted:
        sys.stderr.write("Skipped " + file + " - wrong kmer size\n")
        continue
    if id not in bacteria:
        continue

    thesecounts={}
    with open(os.path.join(dir, file), 'r') as fin:
        for line in fin:
            line = line.strip()
            p = line.split("\t")
            thesecounts[p[0]] = float(p[1])
            if p[0] not in allkmers:
                allkmers[p[0]]=1
    S=sum(thesecounts.values())*1.0
    thesecounts = {x:thesecounts[x]/S for x in thesecounts}
    count[id]=thesecounts

kmers = allkmers.keys()
kmers.sort()

orgs = sorted(bacteria.keys(), key=bacteria.get)
group=0
lastorg=""
with open(outf, 'w') as out:
    out.write("ID\tOrganism\tCode\t" + "\t".join(kmers) + "\n")
    n=1
    for id in count:
        thisorg = bacteria[id]
        if lastorg == thisorg:
            n += 1
        else:
            group += 1
            n=1

        out.write(id + "\t" + bacteria[id] + " (" + str(n) + ")\t" + str(group))
        for kmer in kmers:
            if kmer in count[id]:
                out.write("\t" + str(count[id][kmer]))
            else:
                out.write("\t0")
        out.write("\n")
        
