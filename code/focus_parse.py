'''
Parse the output from batch running focus on all the genomes.

We will take a directory and only process those files that end .focus

'''

import os,sys
import re

try:
    indir = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <directory with focus results>")

taxonomy={}
with open('/home3/redwards/phage/host_analysis/all_host_taxid_taxonomy.txt', 'r') as tin:
    for l in tin:
        p=l.strip().split("\t")
        taxonomy[p[0]]=p[2]

score={}
allmg={}
for f in os.listdir(indir):
    if not f.endswith('.focus'):
        continue
    with open(os.path.sep.join([indir, f])) as fin:
        l = fin.readline()
        if not l.startswith('Input:'):
            sys.stderr.write("ERROR: " + os.path.sep.join([indir, f]) + " does not start with a line that begins with Input:. Is this a modified focus file?\n")
            break
        m=re.findall('fasta/(.*)\.gz', l)
        if m == []:
            sys.stderr.write("ERROR: Can't find a metagenome in " + l + "\n")
            break
        mg=m[0]
        mg=mg.replace('.fna', '')
        mg=mg.replace('.fasta', '')
        allmg[mg]=allmg.get(mg, 1)
        keep=False
        for l in fin:
            if 'Printed the results for the' in l:
                break
            m=re.findall('(\d+)\s+(NC_\d+)\s+(\d+\.\d+)', l)
            if m==[]:
                continue
            if m[0][1] not in score:
                score[m[0][1]]={}
            score[m[0][1]][mg]=m[0][2]


mgs=allmg.keys()
mgs.sort()
print "ID\t" + "\t".join(mgs) + "\tTaxonomy"
for nc in score:
    sys.stdout.write(nc)
    for mg in mgs:
        sys.stdout.write("\t" + score[nc].get(mg, "0"))
    sys.stdout.write("\t" + taxonomy[nc] + "\n")

    
