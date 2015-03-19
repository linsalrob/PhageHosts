'''
Add the scores for every nc/nc pair to text files that have two NC ids
in the first column.
'''


import os 
import sys
from scoring import scoring


try:
    f=sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <file of nc/nc> ")

score = scoring()
wanted = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

first = True
with open(f, 'r') as fin:
    for l in fin:
        l=l.strip()
        if first and not l.startswith('NC'):
            print(l + "\t" + "\t".join(wanted))
            first = False
            continue
        first = False
        p=l.split("\t")
        if not p[0].startswith('NC'):
            sys.stderr.write("Phage " +  p[0] + " does not seem to be an NC id\n")
            continue
        if not p[1].startswith('NC'):
            sys.stderr.write("Bacteria " +  p[1] + " does not seem to be an NC id\n")
            continue
        res = score.score_NC(p[0], p[1])
        output=[]
        for w in wanted:
            if w in res and res[w]:
                output.append("1")
            else:
                output.append("0")
        print(l + "\t" + "\t".join(output))



