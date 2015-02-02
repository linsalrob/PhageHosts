""" 
Count the number of matches between a phage and some hosts for each file
and report the mean, median, and stdev of the number of hits per phage.

You can also do this with perl and bash on the command line:
    perl -lane 'print "$#F\n"' NCIDS_FILE | awk '{s+=$1}END{print "Average:",s/NR}'

"""

import os 
import sys
import rob

try:
    f=sys.argv[01]
except:
    sys.exit(sys.argv[0])


data=[]
with open(f, 'r') as inf:
    for l in inf:
        p=l.strip().split("\t")
        discard = p.pop()
        data.append(len(p))

print("\t".join([os.getcwd(), f, str(rob.mean(data)), str(rob.median(data)), str(rob.stdev(data))]))
