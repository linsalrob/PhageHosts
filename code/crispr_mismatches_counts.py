"""
Generate a list of NC ids of the phages and their hosts for the CRISPR
searches. We accept two different parameters, the number of mismatches
and the BLAST output file. Note that we expect the blast output file
to incldue qlen slen as the last two columns
"""

import rob
import sys
import os
import re


try:
    mm=int(sys.argv[1])
    bf=sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <number of mismatches> <crispr blast output file>. You probably want to use crispr.vs.genomes.blastn-nz\n")


hits={}
with open(bf, 'r') as bin:
    for l in bin:
        p=l.strip().split("\t")

        p[2]=float(p[2])
        for i in range(3, 10):
            p[i]=int(p[i])
        p[10] = float(p[10])
        p[11] = float(p[11])
        p[12] = int(p[12])
        p[13] = int(p[13])

        # mismatches in the sequences is the sum of the difference in
        # length and the mismatches reported by blast
        mismatches = p[4] + (p[13]-p[3])
        if mismatches > mm:
            continue

        m=re.findall('NC_\d+', p[0])
        if m==[]:
            sys.stderr.write("No phage found in " + p[0] + "\n")
            continue
        if len(m) > 1:
            sys.stderr.write("More than one phage found in " + p[0] + "\n")
            continue
        phage = m[0]

        m=re.findall('NC_\d+', p[1])
        if m==[]:
            sys.stderr.write("No host found in " + p[0] + "\n")
            continue
        if len(m) > 1:
            sys.stderr.write("More than one host found in " + p[0] + "\n")
            continue
        host = m[0]

        if phage not in hits:
            hits[phage]={}

        if host not in hits[phage]:
            hits[phage][host]=set()

        hits[phage][host].add(p[2])



for p in hits:
    for h in hits[p]:
        if len(hits[p][h]) > 1:
            print(str(len(hits[p][h])) + " matches from phage " + p + 
                  " to host " + h)
