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
host = phage.phageWithNHosts(5)

try:
    dir  = sys.argv[1]
    outf = sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <directory of kmer counts> <file to write output table to>")

count={}
allkmers={}
organismId={}
for file in os.listdir(dir):
    match = re.findall('NC_\d+', file)
    id = match[0]
    if id not in host:
        sys.stderr.write("Found a sequence with id " + id + " but we don't have enough genomes for it\n")
        continue

    gs = host[id]
    if gs not in count:
        count[gs]=[]
        organismId[gs]=[]

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
    count[gs].append(thesecounts)
    organismId[org].append(id)

kmers = allkmers.keys()
kmers.sort()

group=0
with open(outf, 'w') as out:
    out.write("ID\tOrganism\tCode\t" + "\t".join(kmers) + "\n")
    for org in count:
        n=1
        group += 1
        for i in range(len(count[org])):
            out.write(organismId[org][i] + "\t" + org + " (" + str(n) + ")\t" + str(group));
            n += 1
            for kmer in kmers:
                if kmer in count[org][i]:
                    out.write("\t" + str(count[org][i][kmer]))
                else:
                    out.write("\t0")
            out.write("\n")
        


