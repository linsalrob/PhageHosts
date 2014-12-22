'''Score a blast output file so that each query gets 1/n hits assigned to the database'''

import sys
import os

try:
    bf=sys.argv[1]
    type = sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <blast output file> <all for all hits or best for best hits only>")

# we are going to check for the NC_ids to make sure they are all ones we need
#

tax={}
with open('/home3/redwards/phage/host_analysis/all_host_taxid.txt', 'r') as tin:
    for l in tin:
        p=l.strip().split("\t")
        tax[p[0]]=p[1]

hits={}
n={}
with open(bf, 'r') as bfi:
    for line in bfi:
        p=line.strip().split("\t")
        # a check added for the blastx code. We now assume that the query adn subject are just NC_ ids
        if p[0] not in tax:
            # sys.stderr.write("WARNING " + p[0] + " does not appear to be a valid id. Skipped\n")
            continue
        if p[1] not in tax:
            # sys.stderr.write("Warning (minor): database " + p[1]  + "does not appear to be a valid database sequences.Skipped\n")
            continue
        if p[0] not in hits:
            hits[p[0]]={}
            n[p[0]]=0
        
        hits[p[0]][p[1]]=hits[p[0]].get(p[1], 0) + 1
        n[p[0]]+=1



for i in hits:
    s=sorted(hits[i], key=hits[i].get)
    max=1.0*hits[i][s[0]]/n[i]
    for j in s:
        if type.lower() == "best" and 1.0*hits[i][j]/n[i] != max:
            continue
        print i + "\t" + j + "\t" + str(1.0*hits[i][j]/n[i])

