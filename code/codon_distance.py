'''
Calculate the distance between two codon usages. 
We have two files, the first with just the phages and the second
with their hosts. Then we need to calculate which of the hosts is
closest
'''

import os
import sys
sys.path.append('/home3/redwards/bioinformatics/Modules')
import numpy as np
import scipy


remove_ambiguous = True # do we want ambiguous bases or not
codons = set([
    'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 
    'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 
    'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
    'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 
    'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
    'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 
    'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
    'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'
])

def distance(x, y):
    '''
    Calculate the Euclidean distance between codon usages. An alternate 
    solution would be to use either np.linalg.norm or 
    scipy.spatial but neither of these are working on my system'''
    return np.sqrt(np.sum((x-y)**2))


def remove_ambiguous_bases(header, cds):
    '''
    Remove any codons that contain ambiguous bases.
    '''
    temp=[cds[0]]
    for i in range(1,len(header)):
        if header[i] in codons:
            temp.append(cds[i])

    return temp

try:
    phageF = sys.argv[1]
    bactF = sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <phage file> <hosts file>\n")

cds = {}
header = None
with open(bactF, 'r') as bf:
    for line in bf:
        if line.startswith('Locus'):
            header = line.strip().split("\t")
            for i in range(len(header)):
                header[i] = header[i].strip()
            continue
        line = line.rstrip()
        p = line.split("\t")
        if remove_ambiguous:
            p = remove_ambiguous_bases(header, p)
        cds[p[0]] = np.array([float(x) for x in p[1:len(p)]])


header = None
with open(phageF, 'r') as ph:
    for line in ph:
        if line.startswith('Locus'):
            header = line.strip().split("\t")
            for i in range(len(header)):
                header[i] = header[i].strip()
            continue
        line = line.rstrip()
        p = line.split("\t")
        lowestScore = 1000
        bestHits = []
        if remove_ambiguous:
            p = remove_ambiguous_bases(header, p)
        a1 = np.array([float(x) for x in p[1:len(p)]])
        for c in cds:
            #dist = scipy.spatial.distance.cdist(a1, cds[c])
            #dist = np.linalg.norm(a1-cds[c])
            dist = distance(a1, cds[c])
            if dist < lowestScore:
                lowestScore = dist
                bestHits = [c]
            elif dist == lowestScore:
                bestHits.append(c)
        print p[0]+ "\t" + "\t".join(bestHits)


