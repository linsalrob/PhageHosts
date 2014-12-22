'''Calculate the distance between two codon usages. We will have two files, the first with just the phages and the second with their hosts. Then we need to calculate which of the hosts is closest.

I am worried that ~300 similarities are to NC_014655, and so I will just print out the distances for one phage.'''

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
    phageToPrint = sys.argv[3]
except:
    sys.exit(sys.argv[0] + " <phage file> <hosts file> <phage to print>\n")

cds = {}
with open(bactF, 'r') as bf:
    for line in bf:
        if line.startswith('Locus'):
            continue
        line = line.rstrip()
        p = line.split("\t")
        cds[p[0]] = np.array([float(x) for x in p[1:len(p)]])


with open(phageF, 'r') as ph:
    for line in ph:
        if line.startswith('Locus'):
            continue
        if not line.startswith(phageToPrint):
            continue
        line = line.rstrip()
        p = line.split("\t")
        lowestScore = 1000
        bestHits = []
        a1 = np.array([float(x) for x in p[1:len(p)]])
        for c in cds:
            dist = distance(a1, cds[c])
            print "\t".join([p[0], c, str(dist)])


