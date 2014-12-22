'''
Convert a list with NC_ to a list with NC_ genus species

'''

import sys
sys.path.append('/home3/redwards/bioinformatics/Modules')
import taxon

try:
    f = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <file to convert>")

taxa = taxon.readNodes()
names,blastname,genbankname,synonym = taxon.extendedNames()
divs = taxon.readDivisions()

# read the tax id file
bact={}
with open('/home3/redwards/phage/host_analysis/bacteria_taxid.txt', 'r') as fin:
    for l in fin:
        p=l.strip().split("\t")
        bact[p[0]]=p[1]


with open(f, 'r') as fin:
    for l in fin:
        l=l.strip()
        p=l.split("\t")
        sp=[]
        for a in p:
            if a.startswith("NC_"):
                if a not in bact:
                    sys.stderr.write(a + " was not found in the bacterial file. Skipped\n")
                    continue
                i=bact[a]
                if i not in taxa:
                    sys.stderr.write("ERROR: taxonomy for " + a + " ("+str(i)+") not found\n")
                    continue
                while i != '1' and taxa[i].parent != '1':
                    if taxa[i].rank == 'species':
                        sp.append(names[i].name)
                        i = '1'
                    else:
                        i = taxa[i].parent
        print l + "\t" + "\t".join(sp)





