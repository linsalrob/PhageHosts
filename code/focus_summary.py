'''
some simple stats from the focus analysis

'''

import os,sys

try:
    f=sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <focus parsed output>")

rows=0
colcount={}
nfound=0
nbigone=0
with open(f, 'r') as fin:
    for l in fin:
        p=l.strip().split("\t")
        if p[0] == "ID":
            continue
        cols = len(p)-2
        colcount[cols]=colcount.get(cols, 0) +1
        rows+=1
        found=False
        bo=False
        for s in p[1:-1]:
            if s > 0:
                found=True
            if s >= 1:
                bo=True
        if found:
            nfound += 1
        if bo:
            nbigone+=1

print "Rows: " + str( rows )
print "Column distribution: "
for c in colcount:
    print "\t" + str(c) + ":" + str(colcount[c])

print "Number with any score: " + str(nfound)
print "Number with a score of 1 or more: " + str(nbigone)


