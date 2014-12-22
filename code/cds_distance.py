'''Calculate the distance between two codon usages. We will have two files, the first with just the phages and the second with their hosts. Then we need to calculate which of the hosts is closest'''

import os
import sys
sys.path.append('/home3/redwards/bioinformatics/Modules')
import numpy as np
import scipy


def distance(x, y):
    ''' a better solution would be to use either np.linalg.norm or scipy.spatial but neither of these are working '''
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
        gc[p[0]]=p[-1]

with open(phageF, 'r') as f:
    for line in f:
        if line.startswith('Locus'):
            continue
        p = line.strip().split("\t")
        val = p[-1]
        min = 10000
        best = []
        for h in gc:
            dist=abs(float(gc[h]) - float(val))
            if dist < min:
                min = dist
                best = [h]
            elif dist == min:
                best.append(h)

        print p[0] + "\t" + "\t".join(best)


