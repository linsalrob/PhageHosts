'''Some common functions we need for our phage host prediction work'''
import random
import re
from rob import readFasta

class Phage:
    def __init__(self):
        self.host = {}
        self.hostTaxString={}
        self.hostTax={}
        self.bacteria = {}
        self.fasta = {}

    def phageHost(self):
        '''Get a hash of all phages and the host that phage infects.'''
        if len(self.host) > 0:
            return self.host
        with open('/home3/redwards/phage/host_analysis/phage_host.tsv', 'r') as phin:
            for line in phin:
                line = line.strip()
                [id, gs]=line.split("\t")
                if id in self.host:
                    sys.stderr.write("WARNING: " + id + " was found twice in /home3/redwards/phage/host_analysis/phage_host.tsv\n")
                self.host[id] = gs
        return self.host

    def phageIDs(self):
        '''Return a sorted list of the phage IDs'''
        if len(self.host) == 0:
            self.phageHost()
        ids=self.host.keys()
        ids.sort()
        return ids

    def phageTaxonomyString(self):
        '''Return a hash where the key is the phages we are interested in and the value is the taxonomy string of the host for that phage'''
        phages = self.phageHost()
        with open('/home3/redwards/phage/host_analysis/all_host_taxid_taxonomy.txt', 'r') as fin:
            for line in fin:
                p=line.strip().split("\t")
                if p[0] in phages:
                    self.hostTaxString[p[0]]=p[2]
        return self.hostTaxString

    def phageTaxonomy(self):
        '''Return a hash where the key is the phages we are interested in and the value is the taxonomy id of the host for that phage'''
        phages = self.phageHost()
        with open('/home3/redwards/phage/host_analysis/all_host_taxid_taxonomy.txt', 'r') as fin:
            for line in fin:
                p=line.strip().split("\t")
                if p[0] in phages:
                    self.hostTax[p[0]]=p[1]
        return self.hostTax



    def phageWithNHosts(self, n):
        '''Return a subset of the phages where the host that phage infects has n phages that can infect it. If a host has more than n phages, then the phages are sampled at random.'''
        if len(self.host) == 0:
            self.phageHost()

        self.hostPhages = {}
        for p in self.host:
            self.hostPhages.setdefault(self.host[p], [])
            self.hostPhages[self.host[p]].append(p)
        
        toReturn = []
        for h in self.hostPhages:
            if len(self.hostPhages[h]) < n:
                continue
            if len(self.hostPhages[h]) > n:
                random.shuffle(self.hostPhages[h])
                self.hostPhages[h] = self.hostPhages[h][0:n-1]
            toReturn.extend(self.hostPhages[h])
        return {x:self.host[x] for x in toReturn}
        
    def completeBacteria(self, version=False):
        '''Return a hash of all the complete bacteria that we will use in this study. The key is the ID, and the value is the name of the bacteria. If you want the version number set version to True'''
        if len(self.bacteria) > 0:
            return self.bacteria
        with open('/lustre/usr/data/NCBI/RefSeq/bacteria/complete_genome_ids.txt', 'r') as ids:
            for line in ids:
                line = line.strip()
                source, idstring, name = line.split("\t")
                m = re.findall('ref\|[\w\.]+', idstring)
                if m == None:
                    sys.stderr.write("No refseq id was found in " + idstring + ". Skipped\n")
                    continue
                k=""
                if version:
                    k=m[0].replace('ref|', '')
                else:
                    n=re.findall('ref\|(.*?)\.\d+', m[0])
                    k=n[0]
                self.bacteria[k]=name
        return self.bacteria

    def completeBacteriaIDs(self, version=False):
        '''Get a sorted list of all the bacterial IDs that we use'''
        if len(self.bacteria) == 0:
            self.completeBacteria()
        ids=self.bacteria.keys()
        ids.sort()
        return ids

    def phageSequences(self, fafile='/home3/redwards/phage/host_analysis/phage_with_host.fna'):
        '''Get the DNA sequences of the phages that we are interested in '''
        if len(self.host) == 0:
            self.phageHost()
        
        fa = readFasta(fafile)
        for i in fa:
            m=re.findall('(NC_\d+)', i)
            if m[0] not in self.host:
                continue
            self.fasta[m[0]]=fa[i]
        return self.fasta

    def phageSequenceLengths(self, fafile='/home3/redwards/phage/host_analysis/phage_with_host.fna'):
        '''return a hash with the lengths of the sequences'''
        if len(self.fasta) == 0:
            self.phageSequences(fafile)
        lengths={i:len(self.fasta[i]) for i in self.fasta}
        return lengths
