
import os,sys
# we need to ignore nan's in our values!
from numpy import nanmax
from numpy import nanmin

try:
    df = sys.argv[1]
    type = sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <distance file> <min or max: use min for distance and max for correlation>")

if type not in ["min", "max"]:
    sys.exit(type + " should be one of either 'min' or 'max' depending if you have a distance matrix or correlation scores")

with open(df, 'r') as inf:
    l=inf.readline()
    header = l.strip().split("\t")
    discard = header.pop(0)
    for l in inf:
        p=l.strip().split("\t")
        phage = p.pop(0)
        r=map(float, p)
        if type=='max':
            val=nanmax(r)
        else:
            val=nanmin(r)
        sys.stdout.write(phage)
        for j in xrange(len(r)):
            if r[j] == val:
                sys.stdout.write("\t" + header[j])
        sys.stdout.write("\n")
        
