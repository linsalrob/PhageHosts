'''Convert my tbl files to fasta files of proteins and DNA sequences. Each sequence has the corresponding ID, which is its locus tag, and also has the genome from where it came'''

import sys
import gzip
import os

try:
    dir=sys.argv[1]
except:
    sys.exit(sys.argv[0] + ' <directory with table files>')

prof = open('refseq_proteins.faa', 'w')
dnaf = open('refseq_orfs.faa', 'w')

for f in os.listdir(dir):
    if not f.endswith('genomic.tbl.gz'):
        continue
    fin=gzip.open(os.path.join(dir, f))
    for l in fin:
        l=l.strip()
        p=l.split('\t')
        locus=None
        for tag in p[12].split(';'):
            if tag.startswith('locus:'):
                locus = tag
                locus = locus.replace('locus:', '')
        if locus == None:
            sys.exit('No locus tag found in ' + l)
        location = "_".join([p[5], p[6], p[7]])
        if p[8] == -1:
            location = "_".join([p[5], p[7], p[6]])
        if p[9] and p[10]:
            prof.write('>' + locus + ' [' + p[1]+ '] [' + location + ']\n' + p[9] + "\n")
            dnaf.write('>' + locus + ' [' + p[0]+ '] [' + location + ']\n' + p[10] + "\n")
    fin.close()
prof.close()
dnaf.close()




