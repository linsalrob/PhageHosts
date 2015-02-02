"""A class to hold the data about scoring taxids and NC ids."""


import sys
import taxon

class scoring:
    """
    We use this class to score a pair of IDs. We can handle taxonomy IDs
    and NC_\d+ IDs (i.e. from genbank).
    """

    def __init__(self):
        sys.stderr.write("Parsing the taxonomy files\n")
        self.taxa = taxon.readNodes()
        self.names, self.blastname, self.genbankname, self.synonym = taxon.extendedNames()
        self.divs = taxon.readDivisions()
        sys.stderr.write("Done parsing the taxonomy\n")

        self.wanted = ['species', 'genus', 'family', 'order', 'class',
                       'phylum', 'superkingdom']

        self.taxonomy = {}
        self.nc2tax={}

    def read_host_file(self, taxfile='/home3/redwards/phage/host_analysis/all_host_taxid.txt'):
        """
        Read a host taxonomy file. This file has tuples of NC id and 
        txonomy ID. The important thing here is that the phage taxonomy
        ID is the taxonomy of the host and not the phage!
        """
        with open(taxfile, 'r') as tin:
            for l in tin:
                p=l.strip().split("\t")
                self.nc2tax[p[0]]=p[1]
    
    def score_taxid(self, phage, host):
        """
        Given a pair of taxonomy IDs from the phage and its predicted
        host we return a hash of taxonomic levels with True if the host
        is correct for the phage at that level and False if it is not 
        correct at that level.
        """
        
        results = {w:False for w in self.wanted}

        if phage not in self.taxa:
            sys.stderr.write("SKIPPED Phage: " + str(phage) + " not in the taxonomy\n")
            return results
        if host not in self.taxa:
            sys.stderr.write("SKIPPED host: " + str(host) + " not in the taxonomy\n")
            return results


        if phage not in self.taxonomy:
            self.taxonomy[phage]={}
            i=phage
            while self.taxa[i].parent != '1' and i != '1':
                rank = self.taxa[i].rank
                if rank in self.wanted:
                    self.taxonomy[phage][rank]=i
                i=self.taxa[i].parent

        if host not in self.taxonomy:
            self.taxonomy[host]={}
            i=host
            while self.taxa[i].parent != '1' and i != '1':
                rank = self.taxa[i].rank
                if rank in self.wanted:
                    self.taxonomy[host][rank]=i
                i=self.taxa[i].parent
        
        for w in self.wanted:
            if w in self.taxonomy[phage] and w in self.taxonomy[host]:
                if self.taxonomy[phage][w] == self.taxonomy[host][w]:
                    results[w]=True
        return results

    def score_NC(self, phage, host):
        '''
        Given a pair of NC ids (from GenBank) we convert them to 
        taxonomy IDs and then socre those
        '''
        
        if self.nc2tax == {}:
            self.read_host_file()

        if phage not in self.nc2tax:
            sys.stderr.write("SKIPPED: phage ID: " + phage + "  not in 2tax file\n")
            return self.score_taxid(None, None)
        if host not in self.nc2tax:
            sys.stderr.write("SKIPPED: host ID: " + host + "  not in 2tax file\n")
            return self.score_taxid(None, None)

        return self.score_taxid(self.nc2tax[phage], self.nc2tax[host])


if __name__ == '__main__':
    sys.stderr.write("Initializing ... ")
    s=scoring()

    print("Checking tax id's. These are the same:")
    r = s.score_taxid('1280', '1280')
    print("\n".join([w + "\t" + str(r[w]) for w in r]))
    print("Checking tax id's. These are different")
    r = s.score_taxid('1280', '28901')
    print("\n".join([w + "\t" + str(r[w]) for w in r]))
    print("Checking NC ids. These are the same:")
    r = s.score_NC('NC_021774', 'NC_021775')
    print("\n".join([w + "\t" + str(r[w]) for w in r]))
    print("Checking NC id's. These are different")
    r = s.score_NC('NC_019513', 'NC_021774')
    print("\n".join([w + "\t" + str(r[w]) for w in r]))
