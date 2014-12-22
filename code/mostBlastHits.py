import os
import sys
import re

# read the phage taxonomy file

try:
    bf=sys.argv[1]
except:
    sys.exit(sys.argv[0] + "<blast output file>")

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
        if p[0] not in tax:
            continue
        if p[1] not in tax:
            continue
        if p[0] not in results:
            results[p[0]]={}
        results[p[0]][tax[p[1]]]=results[p[0]].get(tax[p[1]], 0) + 1

sys.stderr.write("We have a total of " + str(len(results)) + " hits\n")
for ph in results:
    k = sorted(results[ph], key=results[ph].get, reverse=True)
    # print ph + "\t" + "\t".join(results[ph].keys())
    # print ph + "\t" + k[0] + "\t" + str(results[ph][k[0]])
    print tax[ph] + "\t" + k[0]
    



