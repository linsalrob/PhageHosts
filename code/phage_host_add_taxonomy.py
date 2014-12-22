import os,sys
import taxon

wanted = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
wanted.reverse()

try:
    inputF = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <input file>")

taxa = taxon.readNodes()
names,blastname,genbankname,synonym = taxon.extendedNames()
divs = taxon.readDivisions()


with open(inputF, 'r') as inf:
    for line in inf:
        if line.startswith('#'):
            print line.strip() + "\tTaxonomy"
            continue
        phage, hosttid = line.strip().split("\t")
        i=hosttid
        taxonomy={}
        while taxa[i].parent != '1' and i != '1':
            rank = taxa[i].rank
            if rank in wanted:
                taxonomy[rank]=names[i].name
            i=taxa[i].parent
        s = "; ".join([taxonomy.get(x, "") for x in wanted])
        print phage + "\t" + hosttid + "\t" + s


