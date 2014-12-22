'''Add the taxonomy id to the refseq list of complete genomes using gi2taxid'''

import os
import sys
import re

# start by reading the gi2taxid file
#

t={}
## add a couple that we don't have for some weird reason
t['284800255']='653938'
t['397678256']='1198627'
with open('/home/db/taxonomy/gi_taxid_nucl.dmp', 'r') as gif:
    for line in gif:
        line = line.strip()
        gi, taxid = line.split("\t")
        t[gi]=taxid



# now read the complete genomes file
with open('complete_genome_ids_taxid.txt', 'w') as out:
    with open('complete_genome_ids.txt', 'r') as cgi:
        for line in cgi:
            line = line.strip()
            p=line.split("\t")
            m = re.findall('gi\|(\d+)\|', p[1])
            if m == None:
                sys.stderr.write("No gi found in " + p[1] + "\n")
                continue
            if m[0] in t:
                out.write(line + "\t" + t[m[0]] + "\n")
            else:
                out.write(line + "\n")
                sys.stderr.write("GI " + m[0] + " not found in /home/db/taxonomy/gi_taxid_nucl.dmp\n")
