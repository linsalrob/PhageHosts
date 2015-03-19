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

lens=phage.phageSequenceLengths()
bctG = set(phage.completeBacteriaIDs())
phgG = set(phage.phageIDs())

for p in phgG:
    count[p]={}

sys.stderr.write("Reading " + f + "\n")
with open(f, 'r') as bin:
    for l in bin:
        p=l.strip().split("\t")
        if p[0] not in phgG:
            continue
        if p[1] not in bctG:
            continue

        if p[1] not in count[p[0]]:
            count[p[0]][p[1]]=[]
            for i in range(lens[p[0]]+1):
                count[p[0]][p[1]].append(0)

        s = int(p[6])
        e = int(p[7])
        if e < s:
            (s,e)=(e,s)
        for i in range(s,e+1):
            count[p[0]][p[1]][i]=1

sys.stderr.write("Summing and printing\n")
for p in count:
    for b in count[p]:
        c = sum(count[p][b])
        print("\t".join([p, b, str(c)]))

