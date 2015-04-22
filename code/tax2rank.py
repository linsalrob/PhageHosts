"""
Convert all the taxonomy IDs to rank so that we are always talking
about the same thing
"""

import sys
import taxon

try:
    rank = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " rank to use ... eg. species or genus")

taxa = taxon.readNodes()



with open("all_host_taxid.txt", 'r') as rin:
    for l in rin:
        p=l.strip().split("\t")
        ori = p[1]
        while taxa[p[1]].rank != rank and p[1] != '1':
            p[1] = taxa[p[1]].parent
        if p[1] == 1:
            sys.exit("Did not find " + rank + " for " + ori)
        print(p[0] + "\t" + p[1])
