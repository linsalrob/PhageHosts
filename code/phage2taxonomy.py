'''Add the taxonomy id to the phage genomes. 

NOTE: WE USE THE TAXONOMY OF THE HOST NOT THE TAXONOMY OF THE PHAGE!!

'''

import os
import sys
import re
sys.path.append('/home3/redwards/bioinformatics/phage_host')
sys.path.append('/home3/redwards/bioinformatics/Modules')
from phage import Phage
import taxon


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
    print "\t".join([id, str(i)])
