import sys
sys.path.append('/home3/redwards/bioinformatics/phage_host')
sys.path.append('/home3/redwards/bioinformatics/Modules')
from phage import Phage
import re
import os
import taxon

''' Code to add all the phage hosts taxonomic heirarchy to the phage host files'''

wanted = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']



manual = {'Acinetobacter genomosp.' : '471', 'Actinobacillus actinomycetemcomitans' : '714', 'alpha proteobacterium' : '34025', 'Bacillus clarkii' : '79879', 'Brevibacterium flavum' : '92706', 'Celeribacter sp.' : '875171', 'Escherichia sp.' : '237777', 'Geobacillus sp.' : '340407', 'Gordonia rubropertincta' : '36822', 'Iodobacter sp.' : '641420', 'Listeria sp.' : '592375', 'Marinomonas sp.' : '127794', 'Methanobacterium thermoautotrophicum' : '145262', 'methicillin-resistant Staphylococcus' : '1280', 'Nitrincola sp.' : '459834', 'Persicivirga sp.' : '859306', 'Salisaeta sp.' : '1392396', 'Sulfitobacter sp.' : '191468'}

phage = Phage()
host = phage.phageHost()

taxa = taxon.readNodes()
names,blastname,genbankname,synonym = taxon.extendedNames()
divs = taxon.readDivisions()


name2id = {names[x].name:x for x in names}
name2id.update({blastname[x].name:x for x in blastname})
name2id.update({genbankname[x].name:x for x in genbankname})
name2id.update({synonym[x].name:x for x in synonym})

for id in host:
    if host[id] in manual:
        i = manual[host[id]]
    else:
        host[id] = host[id].replace(';', '')
        if host[id] not in name2id:
            sys.stderr.write("'" + host[id] + "' : ' xxx '\n")
            continue
        i = name2id[host[id]]

    name = {}
    while taxa[i].parent != '1' and i != '1':
        
        bn=names[i].name
        if i in blastname:
            bn=blastname[i].name

        rank = taxa[i].rank
        if rank in wanted:
            name[rank]=bn
        i=taxa[i].parent
        
    sys.stdout.write(id + "\t" + host[id])
    for w in wanted:
        sys.stdout.write("\t" + name.get(w, ""))
    sys.stdout.write("\n")

