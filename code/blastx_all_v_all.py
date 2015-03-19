'''
Count the number of proteins in common between the phages and the
bacteria based on blastx
'''


import sys
from phage import Phage
phage=Phage()

try:
    f=sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <blast output file converted to NC/NC format. Probably phage.genomes.blastx")

count={}


bctG = phage.completeBacteriaIDs()
phgG = phage.phageIDs()

for p in phgG:
    count[p]={}
    for b in bctG:
        count[p][b]=0

with open(f, 'r') as bin:
    for l in bin:
        p=l.strip().split("\t")
        if p[0] in count and p[1] in count[p[0]]:
            count[p[0]][p[1]] = count[p[0]].get(p[1], 0) + 1

for p in count:
    for b in count[p]:
        print("\t".join([p, b, str(count[p][b])]))





