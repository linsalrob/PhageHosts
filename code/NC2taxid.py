''' Convert a table with two or more columns where the columns are the names of the phages and their hosts to a series of tax ids'''

import sys
import os
import re

try:
    resultsF = sys.argv[1]
except:
    sys.exit( sys.argv[0] + " <results file (tab separated)> <optional turn on verbose output>" )

verbose = False
if len(sys.argv) > 2:
    verbose=True

hostT  = '/home3/redwards/phage/host_analysis/all_host_taxid.txt'

if not os.path.exists(hostT):
    sys.exit(hostT + " does not exist. This is just a two column table of NC id and taxid\n")

taxa={}
with open(hostT, 'r') as hin:
    for line in hin:
        line = line.strip()
        p = line.split("\t")
        taxa[p[0]] = p[1]


with open(resultsF, 'r') as rin:
    for line in rin:
        line = line.strip()
        pieces=line.split("\t")
        for i in range(len(pieces)):
            pieces[i] = pieces[i].strip()
        if pieces[0] not in taxa:
            if verbose:
                sys.stderr.write('>' + pieces[0] + '< not in taxa\n')
            continue
        for i in range(len(pieces)):
            p=pieces[i]
            match = re.findall('NC_\d+', p)
            if match == None or match == []:
                sys.stderr.write("No NC found in " + p + "\n")
                continue
            if match[0] not in taxa:
                sys.stderr.write("Found an NC with no taxonomy: " + match[0] + "\n")
                continue
            if i == 0:
                sys.stdout.write(taxa[match[0]])
            else:
                sys.stdout.write("\t" + taxa[match[0]])
        print



        
