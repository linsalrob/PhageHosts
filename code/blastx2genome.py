'''
Combine the genome information in the file extracted from refseq with the blast output to provide connections to the genomes.
'''

import sys,os,re

try:
    infile=sys.argv[1]
    outfile=sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <blast input file> <blast output file>")

genome={}
with open('/lustre/usr/data/NCBI/RefSeq/bacteria/complete_genome_proteins.faa', 'r') as fin:
    for l in fin:
        m=re.match('>(\S+)\s\[(NC_\d+)', l)
        if m:
            prot, gen=m.groups()
            if prot in genome:
                #sys.stderr.write("WARNING: Duplicate protein sequence: " + prot + "\n")
                genome[prot].append(gen)
            else:
                genome[prot]=[gen]


with open(infile, 'r') as inf:
    with open(outfile, 'w') as out:
        for l in inf:
            p=l.split("\t")
            m=re.findall('(NC_\d+)', p[0])
            if m==None:
                sys.stderr.write("WARNING: No NC number found in " + l)
                continue
            p[0]=m[0]
            if p[1] not in genome:
                sys.stderr.write("WARNING: " + p[1] + " was not found in /lustre/usr/data/NCBI/RefSeq/bacteria/complete_genome_proteins.faa\n")
                continue
            pp = p[1]
            for g in genome[pp]:
                p[1]=g
                out.write("\t".join(p))


