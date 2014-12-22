'''
Generate a list of the best CAIs for each genome.
This takes the data from cai.txt which is basically the output from grep
'''

best={}
all={}
with open('cai.txt', 'r') as fin:
    for line in fin:
        d   = line.strip().split("\t")
        val = float(d[3])
        g   = d[1].replace('genomes/', '')
        p   = d[2].replace('phages/', '')
        if p not in best or best[p] < val:
            best[p] = val
            all[p]  = [g]
        elif best[p] == val:
            all[p].append(g)

for p in all:
    print p + "\t" + "\t".join(all[p])


