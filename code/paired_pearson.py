"""
Calculate pearson correlations between all phages and all hosts in pairs of metagenomes
"""

import os,sys
from scipy.stats.stats import pearsonr
import numpy as np

try:
    bacteriaF = sys.argv[1]
    phageF    = sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <bacterial file> <phage file> Files must have matched columns")

"""
NOTE: We demand that the columns are in the right order because we have pairs
of metagenomes, and we don't a priori know which is paired with which. Therefore
we trust you got them right!
"""

bct={}
with open(bacteriaF, 'r') as bin:
    l=bin.readline()
    bctheaders = l.strip().split("\t")
    for l in bin:
        p=l.strip().split("\t")
        bct[p[0]]=map(float, p[1:])

phg={}
with open(phageF, 'r') as pin:
    l=pin.readline()
    phgheaders = l.strip().split("\t")
    for l in pin:
        p=l.strip().split("\t")
        phg[p[0]]=map(float, p[1:])


# now we just need all pairwise correlations
# provided both are not zero

for p in phg:
    for b in bct:
        assert len(phg[p]) == len(bct[b]), "Phage and bacteria vectors must be the same length"
        newp = []
        newb = []
        for i in range(len(phg[p])):
            if phg[p][i] != 0 and bct[b][i] != 0:
                newp.append(phg[p][i])
                newb.append(bct[b][i])
        if len(newp) < 3:
            # sys.stderr.write("WARNING: Only " + str(len(newp)) + " in common for " +  p + " and " +  b + "\n")
            print("\t".join([p, b, "-"]))
            continue


        pearson, pscore = pearsonr(newp, newb)
        if pearson == np.nan or str(pearson) == "nan":
            # this occurs when we have too few values for a corelation hence the part above
            # sys.stderr.write("Warning NAN for " +  p + " and " +  b + "\n")
            pearson = 0
        print("\t".join([p, b, str(pearson)]))

    

