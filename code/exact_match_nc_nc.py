'''
Generate a list of all pairwise comparisons of the exact matches
'''


import sys
import re
from phage import Phage
phage = Phage()

try:
    f = sys.argv[1]
except:
    sys.exit("Exact match file, probably phage.kmers.bacteria.rc.txt")


bg = phage.completeBacteriaIDs()
pg = phage.phageIDs()


matches={}
for p in pg:
    matches[p]={}
    for b in bg:
        matches[p][b] = 0

with open(f, 'r') as fin:
    for l in fin:
        p=l.strip().split("\t")
        m=re.findall('NC_\d+', l)
        if len(m) != 2:
            #sys.stderr.write("Error parsing two NC ids from " + l)
            continue

        if m[0] not in matches or m[1] not in matches[m[0]]:
            continue

        p[4] = int(p[4])
        if p[4] > matches[m[0]][m[1]]:
            matches[m[0]][m[1]]=p[4]

for p in matches:
    for b in matches[p]:
        print("\t".join([p, b, str(matches[p][b])]))



