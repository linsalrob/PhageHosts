"""
Summarize the scores in the nc/nc files. These are files that have the
phage, the bacteria, some score parameter, and then whether the score is
correct at species, genus, family, order, class, phylum, superkingdom
"""

import sys
import os


usage = '''<file> <max position> <threshold>
The file should be the ncnc.score.tsv file that has the following columns:
    
posistion phage bacteria metric species genus family order class phylum

The max position is the maximum allowable position to include in the score
and we use this value or lower as the max (so use 1 to get only the top
hits).

The threshold is the minimum metric for which we score a successful match.

We use > threshold as the limit (so that you can use 0 as a valid input
to choose non-zero values).
'''

    

try:
    f=sys.argv[1]
    maxposn = int(sys.argv[2])
    threshold = float(sys.argv[3])
except:
    sys.exit(sys.argv[0] + usage)


phage = {}
hostsperphage = {}

taxonomy = ['species', 'genus', 'family', 'order', 'class', 'phylum']
cols = {}
for i in range(len(taxonomy)):
    cols[i+4] = taxonomy[i] 

with open(f, 'r') as fin:
    for l in fin:
        if l.startswith("Phage"):
            continue
        p=l.strip().split("\t")
        if len(p) != 10 and len(p) != 11:
            sys.exit("the line \n" + l + "does not have enough columns. It only has " + str(len(p)) + ". Is this the right file?\n") 
        if p[1] not in phage:
            phage[p[1]]=set()
            hostsperphage[p[1]]=0
        if int(p[0]) > maxposn:
            continue
        if float(p[3]) <= threshold:
            continue
        hostsperphage[p[1]] += 1
        for c in cols:
            if int(p[c]) == 1:
                phage[p[1]].add(cols[c])
            elif int(p[c]) != 0:
                sys.stderr.write("Got neither a 0 nor a 1 for column " + str(c) + " in " + l)


correct = {}
incorrect = {}
for v in cols.values():
    correct[v] = 0
    incorrect[v] = 0

for p in phage:
    for v in cols.values():
        if v in phage[p]:
            correct[v] += 1
        else:
            incorrect[v] += 1

print("Taxonomic level\tCorrect Assignments\tIncorrect Assignments\tTotal Assignments\tPercent Correct")
for t in taxonomy:
    c = "%0.2f" % ((1.0 * correct[t] / (incorrect[t] + correct[t])) * 100)
    print(t + "\t" + str(correct[t]) + "\t" + 
        str(incorrect[t]) + "\t" + str(correct[t] + incorrect[t]) +
        "\t" + str(c))


print("\n")
hpp = hostsperphage.values()
nhpp = "%0.2f" % (1.0 * sum(hpp) / len(hpp))
print("with an average of " + str(nhpp) + " hosts per phage")
