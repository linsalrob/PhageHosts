

import os
import sys
import re
from scoring import scoring

try:
    f=sys.argv[1]
except:
    sys.exit(sys.argv[0] + " file to parse?. Probably phage.kmers.bacteria.txt")

scoringO = scoring()

matches={}
lengths=set()
with open(f, 'r') as fin:
    for l in fin:
        p=l.strip().split("\t")
        m=re.findall('NC_\d+', l)
        if len(m) != 2:
            #sys.stderr.write("Error parsing two NC ids from " + l)
            continue

        p[4] = int(p[4])
        lengths.add(p[4])
        if p[4] > 100:
            p[4] = 100
        if m[0] not in matches:
            matches[m[0]]={}
        if m[1] in matches[m[0]]:
            #sys.stderr.write("Had " +  str(matches[m[0]][m[1]]) + " for " + m[0] + " and " + m[1] + " and now have " + p[4] + "\n")
            if p[4] > matches[m[0]][m[1]]:
                matches[m[0]][m[1]]=p[4]
        else:
            matches[m[0]][m[1]]=p[4]

correct={}
incorrect={}
taxalevels = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
lengths = list(lengths)
lengths.sort()

for l in lengths:
    correct[l]={}
    incorrect[l]={}

allscores={}
for m in matches:
    if m not in allscores:
        allscores[m]={}
    for n in matches[m]:
        if n not in allscores:
            allscores[n]={}
        if n in allscores[m]:
            scores = allscores[m][n]
        else:
            scores = scoringO.score_NC(m, n)
            allscores[m][n] = allscores[n][m] = scores

            
        l=matches[m][n]
        for w in scores:
            if scores[w]:
                correct[l][w] = correct[l].get(w, 0)+1
            else:
                incorrect[l][w] = incorrect[l].get(w, 0)+1

print("Length\t" + "\t".join(taxalevels))
for l in lengths:
    sys.stdout.write(str(l))
    for w in taxalevels:
        if (correct[l].get(w, 0) + incorrect[l].get(w, 0)) == 0:
            percent = 0
        else:
            percent = 100.0 * correct[l].get(w, 0) / (correct[l].get(w, 0) + incorrect[l].get(w, 0))
        sys.stdout.write("\t" + str(percent))
    sys.stdout.write("\n")










