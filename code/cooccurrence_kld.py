"""calculate the KLD as a measure of 'distance' between the bacteria and
the phage. 

    D(P || Q) = -sum(Pi ln(Pi/Qi)) where P = bacteria dist. and Q = phage dist.

For each metagenome we calculate Pi log(Pi/Qi), and then sum that. A
couple of obvious points: Pi and Qi need to be non zero, and KLD is
supposedly less susceptible to normalization.

We assume that the bacteria (Pi) is the real distribution and we are
trying to figure out how far the Phage (Qi) are from the bacteria.

Note that KLD is not symmetric: D(P || Q) != D(Q || P)

"""

import os,sys
from math import log

def calc_kld(p, b):
    assert len(p) == len(b)
    kld = 0
    for i in range(len(b)):
        if p[i] == 0 || b[i] == 0:
            continue
        kld += b[i] * (log(b[i]/p[i]))
    return -kld


try:
    bacteriaF = sys.argv[1]
    phageF    = sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <bacterial file> <phage file>")

bact={}
with open(bacteriaF, 'r') as bin:
    l=bin.readline()
    bactheaders = l.strip().split("\t")
    for l in bin:
        p=l.strip().split("\t")
        taxonomy=p.pop()
        bact[p[0]]=map(float, p[1:])

phage={}
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

        phage[p[0]]=map(float, p[1:])


allbact = bact.keys()
allbact.sort()


## calculate pearson correlations
with open("pearson_correlations.tsv", 'w') as out:
    out.write("Phage\t" + "\t".join(allbact) + "\n")
    for ph in phage:
        out.write(ph)
        for ba in allbact:
            kld = calc_kld(phage[ph], bact[ba])
            out.write("\t" + str(kld))
        out.write("\n")

        
