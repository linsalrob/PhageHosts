import sys
import gzip
import os

try:
    dir=sys.argv[1]
except:
    sys.exit(sys.argv[0] + ' <directory with table files>')

for f in os.listdir(dir):
    if not f.endswith('genomic.tbl.gz'):
        continue
    fin=gzip.open(os.path.join(dir, f))
    for l in fin:
        p=l.split('\t')
        if p[9] and not p[10]:
            print f + "\t" +  p[1]

