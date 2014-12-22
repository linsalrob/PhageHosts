'''Figure out which phages we do not have blastn hits for'''

import os,sys
from phage import Phage
import re

# get a list of all phages
phage=Phage()
phages=phage.phageTaxonomyString()

try:
    blastf=sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <blast file>")

found={}
with open(blastf, 'r') as fin:
    for l in fin:
        p=l.split("\t")
        m=re.findall('(NC_\d+)', p[0])
        found[m[0]]=1

for p in phages:
    if p not in found:
        print "MISSED " + p

