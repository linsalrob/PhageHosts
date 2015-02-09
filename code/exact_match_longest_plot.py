

import os
import sys
import re
from scoring import scoring

try:
    f=sys.argv[1]
except:
    sys.exit(sys.argv[0] + " file to parse?. Probably longest_hits_NCIDs.txt")

scoringO = scoring()

matches={}
lengths=set()
with open(f, 'r') as fin:
    for l in fin:
        p=l.strip().split("\t")

        p[2] = int(p[2])
        if p[2] > 100:
            p[2] = 101
        lengths.add(p[2])
        if p[0] not in matches:
            matches[p[0]]={}
        matches[p[0]][p[1]]=p[2]

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

print("Length\t" + "\t\t\t".join(taxalevels))
print("\t"  + ("\t".join(["Correct", "Incorrect", "", ""])) * len(taxalevels)) 
for l in range(15, 102):
    sys.stdout.write(str(l))
    if l not in correct:
        correct[l]={}
        incorrect[l]={}

    for w in taxalevels:
        if (correct[l].get(w, 0) + incorrect[l].get(w, 0)) == 0:
            percent = 0
        else:
            percent = 100.0 * correct[l].get(w, 0) / (correct[l].get(w, 0) + incorrect[l].get(w, 0))
        #sys.stdout.write("\t" + str(percent))
        sys.stdout.write("\t" + str(correct[l].get(w, 0)) + "\t" + str(incorrect[l].get(w, 0)) + "\t")
    sys.stdout.write("\n")










