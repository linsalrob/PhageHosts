#Raw Data for Phage Host Predictions

This archive contains raw data for phage hosts predictions that includes the raw counts for each of the processes that we have used. For questions about this data, please contact Rob Edwards.

The data is organized in approximately the same manner as the paper by Edwards and Dutilh.

Data files:

- [phage_and_hosts.tsv](phage_and_hosts.tsv): This has the phage NC ID in the first column, the phage name in the second column, and the host genus/species in the third column
- [bacteria_taxid.txt](bacteria_taxid.txt): This file has tuples of bacterial NC ID for the genomes that we used, and taxonomy ID of that genome.
- [all_host_taxid.txt](all_host_taxid.txt): This has an NC ID in the first column and a GenBank taxonomy ID in the second column. If the first column has the ID of a bacteria, then the taxonomy ID is of that bacteria. If the first column has the ID of a phage, then the taxonomy ID is of that phage's host genus and species (as given in phage_and_hosts.tsv). You should generally consider this file to have the "correct" answer.
- [GeneticHomology/blastn/phage_host.blastn](https://edwards.sdsu.edu/PhageHosts/raw_data/GeneticHomology/blastn/phage_host.blastn): This is the raw blastn results of the phage genomes against the complete host genomes
- [GeneticHomology/blastx/phage_host.blastx](https://edwards.sdsu.edu/PhageHosts/raw_data/GeneticHomology/blastx/phage_host.blastx): This is the raw blastx results of the phage genomes against the complete host genomes
- [CRISPR_spacers/pilercr1.06.output.fna](https://edwards.sdsu.edu/PhageHosts/raw_data/CRISPR_spacers/pilercr1.06.output.fna): The output from the [PILE-CR](http://www.biomedcentral.com/1471-2105/8/18) program ran against complete genomes
- [CRISPR_spacers/crispr.vs.genomes.blastn-nz](https://edwards.sdsu.edu/PhageHosts/raw_data/CRISPR_spacers/crispr.vs.genomes.blastn-nz): The raw blastn output of the CRISPR searches against the genomes
- [ExactMatches/phage.bacteria.tsv](https://edwards.sdsu.edu/PhageHosts/raw_data/ExactMatches/phage.bacteria.tsv): Tuples of phage id, bacterial id, start position, sequence, and length of match for all exact matches between phages and bacteria > 15 bp with matches on either strand.
- [OligoNucleotides/GC_CDS/bacteria.gc.tsv](https://edwards.sdsu.edu/PhageHosts/raw_data/OligoNucleotides/GC_CDS/bacteria.gc.tsv): Table of locus, G+C, total bp, and fraction G+C for bacterial genomes
- [OligoNucleotides/GC_CDS/phage.gc.tsv](https://edwards.sdsu.edu/PhageHosts/raw_data/OligoNucleotides/GC_CDS/phage.gc.tsv): Table of locus, G+C, total bp, and fraction G+C for phage genomes
- [OligoNucleotides/CodonUsage/bacteria_codon_counts.tsv](https://edwards.sdsu.edu/PhageHosts/raw_data/OligoNucleotides/CodonUsage/bacteria_codon_counts.tsv): Codon usage counts for bacterial genomes. Note the first row is the sequence for the codon and the order is not guaranteed. You should also ignore any ambiguous codons in the sequences. They are included for completeness, but only make up a tiny fraction of the data.
- [OligoNucleotides/CodonUsage/phage_codon_counts.tsv](https://edwards.sdsu.edu/PhageHosts/raw_data/OligoNucleotides/CodonUsage/phage_codon_counts.tsv): Codon usage counts for phage genomes. Note the first row is the sequence for the codon and the order is not guaranteed. You should also ignore any ambiguous codons in the sequences. They are included for completeness, but only make up a tiny fraction of the data.
- [OligoNucleotides/kmers](https://edwards.sdsu.edu/PhageHosts/raw_data/OligoNucleotides/kmers): contains the raw *k*-mer counts for the bacterial and phage genomes. The kmer files (one per kmer size, obviously), have the NC ID in the first column, the name of the organism in the second column, and a code in the third column, followed by the kmers counts for the headers. The code in the thid column is just a number based on the organism in the second column so that you can sort the data by organisms, but have a unique number for each genus/species.
- [OccurenceProfiles/bacteria_predictions_focus.tsv](https://edwards.sdsu.edu/PhageHosts/raw_data/OccurenceProfiles/bacteria_predictions_focus.tsv): tab separated text of each bacteria in the first column and the abundance of that bacteria in each of the metagenomes as predicted by [FOCUS](http://edwards.sdsu.edu/focus/).
- [OccurenceProfiles/phage_mg_count.tsv](https://edwards.sdsu.edu/PhageHosts/raw_data/OccurenceProfiles/phage_mg_count.tsv): tab separated text of each phage in the first column and the normalized abundance of that phage in each of the metagenomes as calculated from reads that mapped via megablast.

