'''Calculate the distance between two codon usages and report all the 
distances. We will have two files, the first with just the phages
and the second with their hosts. Then we need to calculate which of 
the hosts is closest'''

import os
import sys
sys.path.append('/home3/redwards/bioinformatics/Modules')
import numpy as np
import scipy
from phage import Phage

phage = Phage()
bctG = set(phage.completeBacteriaIDs())
phgG = set(phage.phageIDs())


def distance(x, y):
    ''' a better solution would be to use either np.linalg.norm or 
    scipy.spatial but neither of these are working '''
    return np.sqrt(np.sum((x-y)**2))


try:
    phageF = sys.argv[1]
    bactF = sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <phage file> <hosts file>\n")


gc={}
with open(bactF, 'r') as f:
    for line in f:
        if line.startswith('Locus'):
            continue
        p = line.strip().split("\t")
        if p[0] in bctG:
            gc[p[0]]=p[-1]

sys.stderr.write("Found GC for " + str(len(gc)) + " bacteria\n")

with open(phageF, 'r') as f:
    for line in f:
        if line.startswith('Locus'):
            continue
        p = line.strip().split("\t")
        if p[0] not in phgG:
            continue
        val = p[-1]
        for h in gc:
            dist=abs(float(gc[h]) - float(val))
            print("\t".join([p[0], h, str(dist)]))


