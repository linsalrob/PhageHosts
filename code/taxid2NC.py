'''
Convert a table with one or more columns of taxids to NC ids
'''



import sys
import os
import re

try:
    resultsF = sys.argv[1]
except:
    sys.exit( sys.argv[0] + " <results file (tab separated)>")

hostT  = '/home3/redwards/phage/host_analysis/all_host_taxid.txt'

if not os.path.exists(hostT):
    sys.exit(hostT + " does not exist. This is just a two column table of NC id and taxid\n")

taxa={}
with open(hostT, 'r') as hin:
    for line in hin:
        line = line.strip()
        p = line.split("\t")
        taxa[p[1]] = p[0]


with open(resultsF, 'r') as rin:
    for line in rin:
        line = line.strip()
        pieces=line.split("\t")
        for i in range(len(pieces)):
            pieces[i] = pieces[i].strip()
        for i in range(len(pieces)):
            p=pieces[i]
            if p not in taxa:
                sys.stderr.write("Found a taxonomy with no NC" + p + "\n")
                continue
            if i == 0:
                sys.stdout.write(taxa[match[0]])
            else:
                sys.stdout.write("\t" + taxa[match[0]])
        print



        
