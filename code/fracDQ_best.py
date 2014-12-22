import sys,os

try:
    ff=sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <fracDQ output file>")

fracd={}
bestfracd={}
fracq={}
bestfracq={}
with open(ff, 'r') as fin:
    ## YES: THIS IS BACKWARDS
    for p in fin:
        if p.startswith('#'):
            continue
        l=p.strip().split("\t")
        l[2] = float(l[2])
        l[3] = float(l[3])

        if l[0] not in bestfracq:
            bestfracq[l[0]]=l[2]
            fracq[l[0]]=[]

        if l[0] not in bestfracd:
            bestfracd[l[0]]=l[3]
            fracd[l[0]]=[]

        if l[2] > bestfracq[l[0]]:
            bestfracq[l[0]]=l[2]
            fracq[l[0]]=[]

        if l[3] > bestfracd[l[0]]:
            bestfracd[l[0]]=l[3]
            fracd[l[0]]=[]

        if l[2] == bestfracq[l[0]]:
            fracq[l[0]].append(l[1])

        if l[3] == bestfracd[l[0]]:
            fracd[l[0]].append(l[1])


with open('fracQ.besthits', 'w') as out:
    for p in fracq:
        out.write(p + "\t" + "\t".join(fracq[p]) + "\n")

with open('fracD.besthits', 'w') as out:
    for p in fracd:
        out.write(p + "\t" + "\t".join(fracd[p]) + "\n")
        


