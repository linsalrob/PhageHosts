'''
Split the orf files by genome into a named directory
'''
import rob
import sys
import os
import re

try:
    orfF = sys.argv[1]
    dir  = sys.argv[2]
    type = sys.argv[3]
except:
    sys.exit(sys.argv[0] + " <orfs file> <dir to put them> <type must be one of phage or genome>")

if type not in ['phage', 'genome']:
    sys.exit(type + " is not valid: must be either phage or genome\n")

fa = rob.readFasta(orfF)
if not os.path.exists(dir):
    os.mkdir(dir)

for i in fa:
    if type == 'genome':
        m = re.findall('([\w\.\-]+)\s+\[(\w+)\]\s+\[(.*)\]', i)
        if m == []:
            continue
        gene, genome, locus = m[0]
    else:
        x=i.replace(' COMPLEMENT', '')
        m=x.split(' ')
        gene, genome, locus = m
        with open(os.path.sep.join([dir, genome]), 'a') as out:
            out.write('>' + i + "\n" + fa[i] + "\n")



