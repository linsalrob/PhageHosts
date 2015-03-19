import os,sys
from scipy.stats.stats import pearsonr
import numpy as np

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
sys.stderr.write("Found " + str(len(phage)) + " phages\n")
sys.stderr.write("Found " + str(len(allbact)) + " bacteria\n")


## calculate pearson correlations
with open("pearson_ncnc.tsv", 'w') as out:
    out.write("Phage\tBacteria\tDistance\n")
    for ph in phage:
        for ba in allbact:
            pearson, p = pearsonr(phage[ph], bact[ba])
            if pearson == np.nan or str(pearson)=="nan":
                pearson = 0
            out.write(ph + "\t" + ba + "\t" + str(pearson) + "\n")

        
