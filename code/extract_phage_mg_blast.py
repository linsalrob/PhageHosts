'''
Check a pile of blast output files and a list of genomes and report which genomes were not found
I also print out the results of those we need
'''

import os,sys
import re

try:
    blastdir=sys.argv[1]
    wantedf=sys.argv[2]
    outputdir=sys.argv[3]
except:
    sys.exit(sys.argv[0] + " <blast output directory> <file of ids we want, one per line> <directory to put blast results in>")

want={}
with open(wantedf, 'r') as fin:
    for l in fin:
        want[l.strip()]=1

if not os.path.exists(outputdir):
    os.mkdir(outputdir)

found={}
for f in os.listdir(blastdir):
    if f.endswith('blast') or f.endswith('blast.new') or f.endswith('blastn'):
        with open(os.sep.join([blastdir, f]), 'r') as bin:
            with open(os.sep.join([outputdir, f]), 'w') as out:
                for l in bin:
                    p=l.split('\t')
                    m=re.findall('(NC_\d+)', p[1])
                    if m == []:
                        sys.stderr.write("Warning. Could not find a NC in " + p[1] + "\n")
                    p[1]=m[0]
                    if len(p) < 5:
                        sys.stderr.write("Malformed line in "  + f + "\n" + l + "\n")
                        continue
                    if p[0] in want:
                        out.write(l)
                        found[p[0]]=1
                    elif p[1] in want:
                        out.write(l)
                        found[p[1]]=1


for w in want:
    if w not in found:
        sys.stderr.write("NOT FOUND: " + w + "\n")
