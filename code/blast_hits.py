'''

Report either the best hits, best equal hits or all hits from a  blast search.

FOR THIS VERSION WE ARE USING THE BLAST WITH TAXID ADDED AT THE END OF THE LINE (e.g. from GI to TAX)

Note that we use the HIGHEST bit score here for this calculation. 

We generate a table of [phage NC_ id, host tax id]. We need to restrict these hosts to complete genomes.

'''

import sys
import os
import re

try:
    blastF = sys.argv[1]
    type   = sys.argv[2]
    outF   = sys.argv[3]
except:
    sys.exit(sys.argv[0] + " <blast output file> <type> <output file name to write>\nType can be one of: best (for best hits), equal (for equal best hits), and all (for all hits)")


if type not in ['best', 'equal', 'all']:
    sys.exit(type + " must be one of 'best', 'equal', or 'all'\n")

if not os.path.exists(blastF):
    sys.exit(blastF + " does not exist\n")

if os.path.exists(outF):
    sys.exit(outF + " already exists. NOT OVERWRITING\n")


# read the NC -> taxid file

taxid={}
with open('/home3/redwards/phage/host_analysis/all_host_taxid.txt', 'r') as tin:
    for line in tin:
        t = line.strip().split("\t")
        taxid[t[0]]=t[1]


besthits={}
best={}
with open(blastF, 'r') as fin:
    for line in fin:
        all = line.strip().split('\t')
        gi  = all[0].split('|')
        nc  = re.sub('.\d+$', '', gi[3])
        if nc not in taxid:
            # sys.stderr.write(nc + " was not found in the NC->taxid file all_host_taxid.txt. Skipped\n")
            # this generates a lot of errors now!
            continue
        #nc = taxid[nc]
        bit = float(all[11])
        if nc not in best: 
            best[nc]=bit
            besthits[nc]=[all[14]]
        elif type == "all":
            besthits[nc].append(all[14])
        else:
            if bit > best[nc]:
                best[nc]=bit
                besthits[nc]=[all[14]]
            elif bit == best[nc]:
                besthits[nc].append(all[14])

with open(outF, 'w') as out:
    if type == "all" or type == "equal":
        for nc in besthits:
            out.write(taxid[nc] + "\t" + "\t".join(besthits[nc]) + "\n")
    else:
        sys.stderr.write('Writing because ' + type + "\n")
        for nc in besthits:
            out.write(taxid[nc] + "\t" + besthits[nc][0] + "\n")





