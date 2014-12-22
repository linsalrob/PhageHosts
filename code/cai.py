'''

Calculate the Codon Adaptation Index (CAI) of a genome compared to a reference set of highly expressed genes.

The method here is that described by Xia, X. (2007). An Improved Implementation of Codon Adaptation Index. Evoluationary. Bioinformatics Online 3: 53-58
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2684136

Three rules:
    1. ignore all amino acids only encoded by one codon (M and W in the standard genetic code)
    2. ignore any CAI where a count is 0. Becuase we don't know what the value should be
    3. Split L, R, and S into two groups of codons based on their first two bases

Then we calculate w{i,j} for each codon in a group, which is the frequency of that codon in the group divided by the maximum 

'''

import sys
import os
import math

try:
    refF = sys.argv[1]
    cdsF = sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <refence CDS of highlly expressed genes> <CDS of test genes>")



def readFasta(file, wholeId=False):
    '''Read a fasta file and return a hash. If wholeId is set to false only the first part of the ID (upto the first white space) is returned'''
    try:
        if file.endswith('.gz'):
            f=gzip.open(file, 'rb')
        else:
            f=open(file, 'r')
    except:
        sys.exit("Unable to open file " + file)

    seqs={}
    seq=""
    seqid=""
    for line in f:
        line=line.strip()
        if line.startswith(">"):
            if seqid != "":
                seqs[seqid]=seq
                seq=""
            seqid = line.replace(">", "", 1)
            if not wholeId and seqid.count(" ") > 0:
                seqids = seqid.split(" ")
                seqid = seqids[0]
        else:
            seq = seq + line

    seqs[seqid]=seq
    return seqs

def geneticCode():
    # this is the standard genetic code without any ambiguous bases
    return {
    'TTT' : 'F', 'TCT' : 'S1', 'TAT' : 'Y', 'TGT' : 'C',
    'TTC' : 'F', 'TCC' : 'S1', 'TAC' : 'Y', 'TGC' : 'C',
    'TTA' : 'L0', 'TCA' : 'S1', 'TAA' : '*', 'TGA' : '*',
    'TTG' : 'L0', 'TCG' : 'S1', 'TAG' : '*', 'TGG' : 'W',
    'CTT' : 'L1', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R1',
    'CTC' : 'L1', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R1',
    'CTA' : 'L1', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R1',
    'CTG' : 'L1', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R1',
    'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S0',
    'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S0',
    'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R0',
    'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R0',
    'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G',
    'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G',
    'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G',
    'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G',
    }


def codonGroups():
    # these are the groups that we will use to calculate the codon adaptation index
    # note that M and W have been deleted from this table, and also L, S, and R are 
    # split into two groups each as suggested by Xia.
    return {
        'A': ['GCA', 'GCC', 'GCG', 'GCT'],
        'C': ['TGT', 'TGC'],
        'E': ['GAG', 'GAA'],
        'D': ['GAT', 'GAC'],
        'G': ['GGT', 'GGG', 'GGA', 'GGC'],
        'F': ['TTT', 'TTC'],
        'I': ['ATC', 'ATA', 'ATT'],
        'H': ['CAT', 'CAC'],
        'K': ['AAA', 'AAG'],
        'L0': ['TTA', 'TTG'],
        'L1': ['CTT', 'CTG', 'CTA', 'CTC'],
        'N': ['AAC', 'AAT'],
        'Q': ['CAA', 'CAG'],
        'P': ['CCT', 'CCA', 'CCG', 'CCC'],
        'S0': ['AGC', 'AGT'],
        'S1': ['TCT', 'TCG', 'TCC', 'TCA'],
        'R0': ['AGG', 'AGA'],
        'R1': ['CGA', 'CGG', 'CGT', 'CGC'],
        'T': ['ACA', 'ACT', 'ACG', 'ACC'],
        'V': ['GTA', 'GTC', 'GTG', 'GTT'],
        'Y': ['TAT', 'TAC']
    }

# count the refence highly expressed genes

# get the genetic code and delete the two single use codons. 
# I did it this way so that the above genetic code is intact
gc=geneticCode()
gc.pop('ATG', 0)
gc.pop('TGG', 0)

# the reference sequence
ref = readFasta(refF)
refcount={}

# count the codons in the reference sequence. We do not include ambiguous bases
for i in ref:
    for j in xrange(0, len(ref[i]), 3):
        codon = ref[i][j:j+3]
        if codon == 'ATG' or codon == 'TGG':
            continue
        if codon not in gc:
            sys.stderr.write("In " + refF + " found a codon " + codon + " that is not in the standard genetic code. Ignored\n")
            continue
        refcount[codon] = refcount.get(codon, 0) + 1

# now check to see if any values are zero
ignore = {}
for c in gc:
    if c not in refcount:
        sys.stderr.write("In " + refF + " we have a zero count for " + c + " and so these results will be unreliable!\n")
        ignore[gc[c]] = 1

# calculate w{i,j} for every member of the group
w={}
groups = codonGroups()
for g in groups:
    if g in ignore:
        for c in groups[g]:
            sys.stderr.write("Ignoring " + c + " (" + g + ")\n")
            w[c]=0
    else:
        vals=[]
        for c in groups[g]:
            vals.append(refcount[c])
        maxv = max(vals)
        for c in groups[g]:
            w[c] = 1.0 * refcount[c]/maxv

# now we read the ORFs that we are trying to compare and count those
cds = readFasta(cdsF)
cdscount={}
for i in cds:
    for j in xrange(0, len(cds[i]), 3):
        codon = cds[i][j:j+3]
        if codon == 'ATG' or codon == 'TGG':
            continue
        if codon not in gc:
            sys.stderr.write("In " + cdsF + " found a codon " + codon + " that is not in the standard genetic code. Ignored\n")
            continue
        cdscount[codon] = cdscount.get(codon, 0) + 1


cai = 0
for c in gc:
    if c not in w or c not in cdscount or w[c] == 0:
        # these will add zero to the cai
        continue
    cai += (cdscount[c] * math.log(w[c]))/cdscount[c]

print "CAI\t" + refF + "\t" + cdsF + "\t" + str(math.exp(cai))



        







