'''Count the number of hits to each of the phages in the different metagenomes'''

import sys,os,re
from phage import Phage
import gzip

# get a list of all phages
phage=Phage()
phages=phage.phageTaxonomyString()
count={}
for p in phages:
    count[p] = {}

# read the sizes of all the metagenomes. This will be slow but oh well
nseqs={}
with open('/home3/redwards/phage/host_analysis/cooccurence/blast/metagenome_statistics.tsv', 'r') as min:
    for l in min:
        if l.startswith('ID'):
            continue
        p=l.strip().split("\t")
        nseqs[p[0]]=int(p[1])

if False:
    for f in os.listdir('/home3/katelyn/METAGENOMES/MGRAST/fasta'):
        name=f.replace('.fasta.gz', '')
        name=name.replace('.fna.gz', '')
        nseqs[name]=0
        fin=gzip.open(os.path.join('/home3/katelyn/METAGENOMES/MGRAST/fasta', f))
        for line in fin:
            if line.startswith('>'):
                nseqs[name] +=1
        sys.stderr.write("Read " + f + " and saving as " + name + " had " + str(nseqs[name]) + " sequences\n")
        fin.close()

ignoredphage={}
nseqsk=nseqs.keys()
for m in nseqsk:
    if not os.path.exists(os.path.join('/home3/redwards/phage/host_analysis/cooccurence/blast/combined_blast', m + ".blast")):
        sys.stderr.write("No blast results found for metagenome: " + m + "\n")
        nseqs.pop(m)
        continue
    
    for p in count:
        count[p][m] = 0
    with open(os.path.join('/home3/redwards/phage/host_analysis/cooccurence/blast/combined_blast', m + ".blast"), 'r') as fin:
        for line in fin:
            p=line.split("\t")
            phageId=p[1]
            if not phageId.startswith('NC'):
                matches=re.findall('(NC_\d+)', phageId)
                if matches == []:
                    sys.stderr.write("WARNING: No NC ID found in " + phageId + "\n")
                    continue
                phageId = matches[0]
            if phageId not in count:
                #sys.stderr.write("Found phage " + p[1] + " that we don't want. Skipped\n")
                ignoredphage[phageId]=1
                continue
            count[phageId][m] += 1

sys.stderr.write("Deleted a total of " + str(len(ignoredphage)) + " phages\n")
for p in ignoredphage:
    sys.stderr.write(p + "\n")

mgs=nseqs.keys()
mgs.sort()

with open('raw_phage_mg_counts.tsv', 'w') as raw:
    with open('normalized_phage_mg_counts.tsv', 'w') as norf:
        norf.write("ID\t" + "\t".join(mgs) + "\tTaxonomy\n")
        raw.write("ID\t" + "\t".join(mgs) + "\tTaxonomy\n")
        for p in phages:
            raw.write(p)
            norf.write(p)
            for m in mgs:
                raw.write("\t" + str(count[p][m]))
                norf.write("\t" + str(1.0*count[p][m]/nseqs[m]))
            raw.write("\t" + phages[p] + "\n")
            norf.write("\t" + phages[p] + "\n")

