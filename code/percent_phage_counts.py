"""
The phage counts per metagenome are normalized based on the number of 
reads that hit. I want to scale that to a percent, so that it matches
the bacterial data. If there was a single phage present it would get 
100% of the reads
"""


import os
import sys


try:
    inf = sys.argv[1]
    ouf = sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <phage abundance file> (probably normalized_phage_mg_counts.tsv) <output file>")


header = None
data = []
total = []
with open(inf, 'r') as fin:
    header = fin.readline()
    h = header.split("\t")
    for i in range(len(h)):
        total.append(0)
    for l in fin:
        p = l.strip().split("\t")
        for i in range(1, len(p)-1):
            total[i] += float(p[i])
        data.append(p)

with open(ouf, 'w') as out:
    out.write(header)
    for l in data:
        out.write(l[0])
        for i in range(1, len(l)-1):
            if total[i] == 0:
                # this metagenome has no phage hits!
                out.write("\t0")
            else:
                out.write("\t" + str(1.0 * float(l[i])/total[i]  * 100))
        out.write("\t" + l[-1] + "\n")

