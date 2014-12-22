import os
import sys
import re

# read the phage taxonomy file

try:
    bf=sys.argv[1]
except:
    sys.exit(sys.argv[0] + "<parsed blast output file>")

phage={}
with open('/home3/redwards/phage/host_analysis/phage_taxid.txt', 'r') as p:
    for line in p:
        a=line.strip().split("\t")
        phage[a[0]]=a[1]


# read the crispr id to taxonomy file
crispr={}
with open('/home/db/CRISPR/crispr.u-psud.fr/id.taxa', 'r') as p:
    for line in p:
        a=line.strip().split("\t")
        crispr[a[0]]=a[1]


# read the blast output file

results={}
with open(bf, 'r') as p:
    for line in p:
        parts=line.strip().split("\t")
        # convert phage to taxonomy
        m=re.findall('(NC_\d+)', parts[0])
        ph = phage[m[0]]
        if ph not in results:
            results[ph]={}
        # convert the crispr to taxonomy
        cr = crispr[parts[1]]
        results[ph][cr]=1

for ph in results:
    print ph + "\t" + "\t".join(results[ph].keys())




