''' 
Score the assignments of phage to host.

The input file should have the phage taxonomy id in the first coulmn, and then the host(s) taxonomy id in the subsequent columns. 

(See the code NC to taxid to create these data from a tab separated list of files. This is based on previously calculated taxa information.)

We will report the accuracy at each of the phylogenetic levels shown in the wanted list

'''

import sys
sys.path.append('/home3/redwards/bioinformatics/phage_host')
sys.path.append('/home3/redwards/bioinformatics/Modules')
from phage import Phage
import re
import os
import taxon

try:
    inputF = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <file of taxon ids separated with tabs (phage in first column)>")


wanted = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
match={w:0 for w in wanted}
mismatch={w:0 for w in wanted}

taxa = taxon.readNodes()
names,blastname,genbankname,synonym = taxon.extendedNames()
divs = taxon.readDivisions()


with open(inputF, 'r') as inf:
    for line in inf:
        line = line.strip()
        hosts = line.split('\t')
        phage = hosts.pop(0)
        taxonomy={}
        taxonomy[phage]={}
        for h in hosts:
            taxonomy[h]={}

        if phage not in taxa:
            sys.stderr.write("SKIPPED: Phage " + str(phage) + " not found in the taxonomy\n")
            continue
        i=phage
        while taxa[i].parent != '1' and i != '1':
            rank = taxa[i].rank
            if rank in wanted:
                taxonomy[phage][rank]=i
            i=taxa[i].parent

        best=0
        bestHost=None
        for h in hosts:
            count=0
            i=h
            if i not in taxa:
                sys.stderr.write("ERROR: " + str(i) + " is not in the taxonomy file\n")
                continue
            while taxa[i].parent != '1' and i != '1':
                rank = taxa[i].rank
                if rank in wanted:
                    taxonomy[h][rank]=i
                    if rank in taxonomy[phage] and i == taxonomy[phage][rank]:
                        count += 1
                i=taxa[i].parent
            if count > best:
                best = count
                bestHost = h

        if bestHost == None:
            # none of the hosts had a single similarity to the phage
            for w in wanted:
                mismatch[w] += 1
        else:
            # after a discussion with Bas we decided to not score non-existent taxonomies
            for w in wanted:
                if w not in taxonomy[phage] and w not in taxonomy[bestHost]:
                    #match[w] +=1
                    continue
                elif w not in taxonomy[phage]:
                    # sys.stderr.write("No " + w + " in taxonomy for phage: " + phage + "\n")
                    # mismatch[w] +=1
                    continue
                elif w not in taxonomy[bestHost]:
                    # sys.stderr.write("No " + w + " in taxonomy for host: " + bestHost + "\n")
                    # mismatch[w] +=1
                    continue
                elif taxonomy[phage][w] == taxonomy[bestHost][w]:
                    match[w] +=1
                else:
                    mismatch[w] +=1

print "Taxonomic level\tIncorrect assignments\tCorrect assignments\tTotal Assignments\tPercent Correct"
for w in wanted:
    t = mismatch[w] + match[w]
    f = 1.0 * int((1.0 * match[w]/t)*10000)/100

    print "\t".join([w, str(mismatch[w]), str(match[w]), str(t), str(f)])





