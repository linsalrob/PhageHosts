import os
import sys
import re

# read the phage taxonomy file

try:
    bf=sys.argv[1]
except:
    sys.exit(sys.argv[0] + "<parsed blast output file>")

# read the taxonomy file
tax={}
with open('/home3/redwards/phage/host_analysis/all_host_taxid.txt', 'r') as tin:
    for l in tin:
        p=l.strip().split("\t")
        tax[p[0]]=p[1]

# read the blast output file

results={}
with open(bf, 'r') as p:
    for l in p:
        p=l.strip().split("\t")
        if tax[p[0]] not in results:
            results[tax[p[0]]]={}
        results[tax[p[0]]][tax[p[1]]]=1

for ph in results:
    print ph + "\t" + "\t".join(results[ph].keys())




