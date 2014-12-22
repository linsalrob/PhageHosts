PhageHosts
==========

The code used for identification of the hosts of different phages. This is the complete code base used in Edwards, McNair, Wellington-Oguri, and Dutilh, manuscript in preparation. 

Almost all of the code is written in python and should work with version 2.7. Some parts of the code require matplotlib, numpy, and/or scipy. Some parts of the code were written in Perl 5. The code was written on CentOS 6 machines. Parts of the code include reference to running things on our cluster where we use SGE as the job scheduler. You can run those parts without the cluster, but it will probably be slower!


0. Datasets
===========

The phage genomes were downloaded from genbank and refseq and parsed to get the information in the "host" field. 

Use python ~/bioinformatics/ncbi/combine_gbff_fna.py viral.1.genomic.gbff viral.1.1.genomic.fna > viral.1.cds.fna
to convert from GBFF to FNA format for all of the open reading frames

There are 1046 phages in genbank that have a host annotation:
perl -ne '@a=split /\t/; print if ($a[2])' genbank.txt | grep YES$ | grep -i 'complete genome' | wc -l

There are 971 phages in refseq that have a host annotation:
PERL -ne '@a=split /\t/; print if ($a[2])' refseq.txt | grep YES$ | grep -i 'complete genome' | wc -l

We use the refseq phages because we are going to get the refseq genomes.
use perl get_viral_dna.pl  to split the viruses into either phage or eukaryotic viruses. Obv we need phage for this!

Make a table of the heirarchy of the phage hosts. We are just going to keep these ranks: species, genus, family, order, class, phylum, superkingdom

I wrote code to automatically get the taxonomy for most of the hosts, but there were a few that I couldn't map, so I added those manually. To whit:
	'Acinetobacter genomosp.' : '471',
	'Actinobacillus actinomycetemcomitans' : '714',
	'alpha proteobacterium' : '34025',
	'Bacillus clarkii' : '79879',
	'Brevibacterium flavum' : '92706',
	'Celeribacter sp.' : '875171',
	'Escherichia sp.' : '237777',
	'Geobacillus sp.' : '340407',
	'Gordonia rubropertincta' : '36822',
	'Iodobacter sp.' : '641420',
	'Listeria sp.' : '592375',
	'Marinomonas sp.' : '127794',
	'Methanobacterium thermoautotrophicum' : '145262',
	'methicillin-resistant Staphylococcus' : '1280',
	'Nitrincola sp.' : '459834',
	'Persicivirga sp.' : '859306',
	'Salisaeta sp.' : '1392396',
	'Sulfitobacter sp.' : '191468'
