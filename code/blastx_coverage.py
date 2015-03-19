'''
Calculate the coverage of each genome pair from the blastx results

Starting with the blastx converted to NC/NC ids, we want to calculate
the coverage at each position in every genome.

We will consider all the genomes with the most number of bases in the
phage as the top genomes

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

sys.stderr.write("Found " + str(len(count)) + ' matches\n')
for p in count:
    tot=0
    genomes=[]
    for b in count[p]:
        c = sum(count[p][b])
        if c > tot:
            tot = c
            genomes = [b]
        elif c == tot:
            genomes.append(b)
    print(p + "\t" + "\t".join(genomes))
sys.stderr.write("Done")
