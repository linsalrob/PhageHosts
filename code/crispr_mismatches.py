"""Count the crispr hits from blastn and allow for a certain number of
mismatches. Report the mean, median, and stdev hosts per phage"""

import rob
import sys
import os
import re

# allowable errors == sum(gaps + mismatches)
errors = 1000
# gaps = # gaps
gaps = 0
# mismatches = # mismatches
mismatches = 1000

try:
    bf=sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <crispr blast output file>. You probably want to use crispr.vs.genomes.blastn-nz\n")

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

        # p[4] is gap, p[5] is mismatches
        if p[4] + p[5] > errors:
            continue

        if p[4] > mismatches:
            continue

        if p[5] > gaps:
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
            hits[phage]=set()
        hits[phage].add(host)

results=[]
for p in hits:
    results.append(len(hits[p]))

estr = " no limit on combined gaps and mismatches, "
if errors < 1000:
    estr = " a maximum of " + str(errors) + " combined gaps and mismatches, "
mstr = " no limit on mismatches, and "
if mismatches < 1000:
    mstr = " a maximum of " + str(mismatches) + " mismatches, and "
gstr = " no limit on gaps, "
if gaps < 1000:
    gstr = " a maximum of " + str(gaps) + " gaps, "


print("With" + estr + mstr + gstr + "we predicted a hosts for " + str(len(hits)) + " phages and there were a mean of " +
      str(rob.mean(results)) + " hosts per phage, a median of " + str(rob.median(results)) + 
      " hosts per phage, and a stdev of " + str(rob.stdev(results)) + " hosts per phage\n")





