import os,sys
from scipy.spatial.distance import cdist
import numpy as np

try:
    bacteriaF = sys.argv[1]
    phageF    = sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <bacterial file> <phage file>")

bactids=[]
data=[]
with open(bacteriaF, 'r') as bin:
    l=bin.readline()
    bactheaders = l.strip().split("\t")
    for l in bin:
        p=l.strip().split("\t")
        taxonomy=p.pop()
        bact=map(float, p[1:])
        if len(bact) - bact.count(0) > 2:
            data.append(bact)
            bactids.append(p[0])

bactdata = np.array(data)

data=[]
phageids=[]
with open(phageF, 'r') as bin:
    l=bin.readline()
    phageheaders = l.strip().split("\t")
    # check that the columns are in the same order (hopefully sorted?)
    reorderCols = False
    for i in xrange(len(phageheaders)):
        if phageheaders[i] != bactheaders[i]:
            reorderCols = True
    colorder=[]
    if reorderCols:
        for i in xrange(len(bactheaders)):
            if bactheaders[i] not in phageheaders:
                sys.exit('FATAL column ' + bactheaders[i] + ' was not found in the phages')
            colorder.append(phageheaders.index(bactheaders[i]))
            
    
    for l in bin:
        p=l.strip().split("\t")
        taxonomy=p.pop()
        
        if reorderCols:
            temp=[]
            for i in xrange(colorder):
                temp[i]=p[colorder[i]]
            p=temp

        phage=map(float, p[1:])
        if len(phage) - phage.count(0) > 2:
            data.append(phage)
            phageids.append(p[0])

phagedata = np.array(data)

allbact=[]
for i in xrange(len(bactids)):
    allbact.append(bactids[i])

for type in ["euclidean", "braycurtis", "cityblock", "hamming", "jaccard"]:
    with open(type + "_correlations.tsv", 'w') as out:
        out.write("Phage\t" + "\t".join(allbact) + "\n")
        distance = cdist(bactdata, phagedata, metric=type)
        for i in xrange(len(phageids)):
            out.write(phageids[i])
            for j in xrange(len(bactids)):
                out.write("\t" + str(distance[j][i]))
            out.write("\n")

