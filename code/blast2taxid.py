'''Add the tax id to a blast output file. This will place the tax id on the end of the line, so the rest of the columns are the same'''

import sys
sys.path.append('/home3/redwards/bioinformatics/Modules')
import taxon


try:
    f=sys.argv[1]
    type=sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <blast output file> <type prot or nucl>")

gi2t=taxon.readGiTaxId(type)

with open(f, 'r') as fin:
    for line in fin:
        line = line.strip()
        p  = line.split("\t")
        p2 = p[1].split("|")
        gi = None
        #print p2
        for i in xrange(len(p2)):
            if p2[i] == "gi":
                gi=p2[i+1]
        if not gi:
            sys.stderr.write("No gi found in " + line + ". Skipped\n")
            continue
        if gi not in gi2t:
            sys.stderr.write("gi " + gi + " not found in table. Skipped\n")
            continue
        print line + "\t" + gi2t[gi]
