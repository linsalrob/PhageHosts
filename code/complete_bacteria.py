import sys
import os

want={}
found={}
with open('/home3/redwards/phage/host_analysis/bacteria_taxid.txt', 'r') as f:
    for line in f:
        p=line.strip().split("\t")
        want[p[0]]=1

for f in os.listdir('.'):
    if not f.startswith('bacteria'):
        continue
    with open(f, 'r') as fin:
        for line in fin:
            p=line.split("\t")
            p[0] = p[0].strip()
            if p[0] in want:
                sys.stdout.write(line)
                found[p[0]]=1

for w in want:
    if w not in found:
        sys.stderr.write("No data for " + w + "\n")
