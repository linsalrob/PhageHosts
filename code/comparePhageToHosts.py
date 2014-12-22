'''Check our list of phages and hosts and make sure that we match at genus species, and higher taxonomic levels'''

import sys,os
import taxon

wanted = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
match={w:0 for w in wanted}
mismatch={w:0 for w in wanted}
exactmatch=0

taxa = taxon.readNodes()
names,blastname,genbankname,synonym = taxon.extendedNames()
divs = taxon.readDivisions()

bac={}
tree={r:[] for r in wanted}
with open('bacteria_taxid.txt', 'r') as inf:
    for l in inf:
        p=l.strip().split("\t")
        bac[p[1]]=1
        i=p[1]
        while taxa[i].parent != '1' and i != '1':
            rank = taxa[i].rank
            if rank in wanted:
                tree[rank].append(i)
            i=taxa[i].parent

n=0
with open('phage_host_taxid_allpresent.txt', 'w') as out:
    with open('phage_host_taxid.txt', 'r') as inf:
        for l in inf:
            if l.startswith('#'):
                continue
            n+=1
            p=l.strip().split("\t")
            if p[1] in bac:
                exactmatch+=1

            i=p[1]
            while taxa[i].parent != '1' and i != '1':
                rank = taxa[i].rank
                if rank in wanted:
                    if i in tree[rank]:
                        match[rank]+=1
                        if rank == 'species':
                            out.write(l)
                    else:
                        mismatch[rank]+=1
                        if rank == 'superkingdom':
                            sys.stderr.write("No superkingdom for " + p[1] + " a " + names[i].name + " virus\n")
                        if rank == 'species':
                            sys.stderr.write("No species for " + p[1] + " a " + names[i].name + " virus\n")

                i=taxa[i].parent



print "Out of " + str(n) + " phages, there were " + str(exactmatch) + " exact matches to the host"
for r in wanted:
    print r + " match: " + str(match[r]) + " mismatch: " + str(mismatch[r])
