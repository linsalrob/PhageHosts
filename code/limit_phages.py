'''
Limit out output to just the interesting 820 phages
'''

import os
import sys

want=set()
with open('/home3/redwards/phage/host_analysis/phage_host.tsv', 'r') as pin:
    for l in pin:
        p=l.strip().split("\t")
        want.add(p[0])

try:
    f = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <file with phages in the first column>")

with open(f, 'r') as pin:
    for l in pin:
        p=l.split('\t')
        if p[0] in want:
            print(l.strip())
