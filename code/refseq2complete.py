import sys,os,re

try:
    bacf=sys.argv[1]
    refs=sys.argv[2]
    outf=sys.argv[3]
except:
    sys.exit(sys.argv[0] + " <list of bacteria> <refseq proteins> <output file>")

want={}
with open(bacf, 'r') as rin:
    for l in rin:
        p=l.split("\t")
        want[p[0]]=1
p=False
with open(refs, 'r') as rin:
    with open(outf, 'w') as out:
        for l in rin:
            if l.startswith('>'):
                p=False
                m=re.findall('\[(NC_\d+)[\.\]]', l)
                if m==[]:
                    sys.stderr.write("No NC in " + l)
                elif m[0] in want:
                    p=True
            if p:
                out.write(l)

