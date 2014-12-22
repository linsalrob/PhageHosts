
# choose the backend so we don't try and use X11
import matplotlib
matplotlib.use("AGG")


import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.lda import LDA

import os
import sys
import math


try:
    infile = sys.argv[1]
    outputfig = sys.argv[2]
except:
    sys.exit(sys.argv[0] + "<input file> <figure to create>\n")


def oldcolors(n):
    '''Generate some pretty colors. This is taken from http://stackoverflow.com/questions/10731147/evenly-distributed-color-range-depending-on-a-count'''
    phi=0
    colors=[]
    red=[]
    grn=[]
    blu=[]
    for i in range(n):
        U=math.cos(phi)
        V=math.sin(phi)
        Y=0.5

        red.append(Y+V/0.88)
        grn.append(Y-0.38*U-0.58*V)
        blu.append(Y+U/0.49)
        phi += 360/numToPlot

    # now we normalize and scale to 255 for RGB
    maxR = max(red); minR = min(red); scaleR = 255/(maxR-minR);
    maxG = max(grn); minG = min(grn); scaleG = 255/(maxG-minG);
    maxB = max(blu); minB = min(blu); scaleB = 255/(maxB-minB);
    
    for i in range(n):
        R = int(scaleR * (red[i]-minR))
        G = int(scaleG * (grn[i]-minG))
        B = int(scaleB * (blu[i]-minB))
        colors.append([R,G,B, 1])
    return colors

def colors(n):
    '''An alternative way of choosing colors'''
    colors = cm.rainbow(np.linspace(0, 1, n))
    return colors
        

data=[]
allkmerlist=[]
code=[]
names=[]
ids=[]
with open(infile, 'r') as tm:
    line = tm.readline()
    line = line.strip()
    allkmerlist=line.split("\t")
    allkmerlist = allkmerlist[1:]
    
    for line in tm:
        line = line.strip()
        pieces = line.split("\t")
        ids.append(pieces[0])
        names.append(pieces[1])
        code.append(int(pieces[2]))
        data.append([float(x) for x in pieces[3:]])


print str(data)        
X=np.array(data)
pca = PCA()
X_r = pca.fit(X).transform(X)

allcodes = np.array(range(max(code)+1))
# what are the 10 most abundant organisms
numToPlot=20
orgcount={}
for c in code:
    if c in orgcount:
        orgcount[c]+=1
    else:
        orgcount[c]=1

allcodes = sorted(allcodes,key=orgcount.get)
allcodes.reverse()
allcodes = allcodes[:numToPlot]
colrs = colors(numToPlot)


print('explained variance ratio (all components): %s'
      % str(pca.explained_variance_ratio_))

plt.figure()


for cl, i, name in zip(colrs, allcodes, names):
    plt.scatter(X_r[code == i, 0], X_r[code == i, 1], c=cl, label=name, )
#plt.legend()
plt.title('PCA of K-mer dataset ' + infile)
#plt.show()
plt.savefig(outputfig)
