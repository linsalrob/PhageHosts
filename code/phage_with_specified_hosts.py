'''

Print a list of phage genomes that have a set of host genus/species

'''

import sys
try:
    hF = sys.argv[1]
    pF = sys.argv[2]
    cF = sys.argv[3]
except:
    sys.exit(sys.argv[0] + " <host file with id and genus_species> <phage.host file> < phage cds file>")

wanted={}
with open(hF, 'r') as hfi:
    for l in hfi:
        p=l.strip().split("\t")
        wanted[p[2]]=1

phage={}
with open('phages_wanted.txt', 'w') as out:
    with open(pF, 'r') as pfi:
        for l in pfi:
            p=l.strip().split("\t")
            if p[1] in wanted:
                phage[p[0]]=1
                out.write(l)

with open('phages_wanted_cds.dna', 'w') as out:
    with open(cF, 'r') as cfi:
        p=0
        for l in cfi:
            if l.startswith('>'):
                p=0
                parts = l.split(" ")
                if parts[1] in phage:
                    p=1
            if p > 0:
                out.write(l)

