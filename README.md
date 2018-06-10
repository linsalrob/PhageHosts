# PhageHosts

The code used for identification of the hosts of different phages. This is the complete code base used in * Robert A. Edwards, Katelyn McNair, Karoline Faust, Jeroen Raes, and Bas E. Dutilh (2015) Computational approaches to predict bacteriophage–host relationships. FEMS Microbiology Reviews [doi: 10.1093/femsre/fuv048](http://dx.doi.org/10.1093/femsre/fuv048) ∗ 

Almost all of the code is written in python and should work with version 2.7. Some parts of the code require matplotlib, numpy, and/or scipy. Some parts of the code were written in Perl 5 and should run with the standard Perl libraries. The code was written on CentOS 6 machines and should run *out of the box* on those machines, but it will also run on Linux, MacOSX, or Windows. Parts of the code include reference to running things on our cluster where we use SGE as the job scheduler. You can run those parts without the cluster, but it will probably be slower!

The data/ directory includes some of the data sets that we used. We have not uploaded all genome sequences: you can get those from [RefSeq](http://www.ncbi.nlm.nih.gov/refseq/).

# Websites

See also [Edwards Lab Website](http://edwards.sdsu.edu/PhageHosts/)

# Datasets

## Phages

The phage genomes were downloaded from RefSeq and parsed to get the information in the "host" field. These data were converted to coding and protein sequences:

For example, to convert from GBFF to FNA format for all of the open reading frames:

```
    python PhageHosts/code/combine_gbff_fna.py viral.1.genomic.gbff viral.1.1.genomic.fna > viral.1.cds.fna
```

There are 971 phages in RefSeq that have a host annotation:

```
    perl -ne '@a=split /\t/; print if ($a[2])' refseq.txt | grep YES$ | grep -i 'complete genome' | wc -l
```

use `perl PhageHosts/code/get_viral_dna.pl`  to split the viruses into either phage or eukaryotic viruses. For this work we are just going to use the Phage datasets. It is left up to the reader to try some of these challenges on Eukaryotic viral data sets.


### Extracting the DNA and protein sequences


To get the phage coding sequences, we pull them out of the fasta file, and then translate them

```
    for i in $(cut -f 1 phage_with_host.tsv); do 
        grep -A 1 $i 2014_07_21/refseq/viral.1.cds.fna; done > data/phage_with_host.cds.fna
        perl PhageHosts/code/translate.pl phage_with_host.cds.fna > data/phage_with_host.cds.faa
    done
```

## Bacteria

All bacterial sequences were downloaded from RefSeq:

```
    ncftpget ftp://ftp.ncbi.nih.gov/refseq/release/bacteria/bacteria.*.fna.gz
```

These were extracted and a single bacterial database was made for blastn searches

```
    cat bacteria.*fna > bacteria.genomic.fna
```

We extract the NC_ ids from the complete bacteria, above:

```
    grep \> bacteria.genomic.fna > ids.txt
    perl -i -npe 's/\:/\t/;s/^(.*)\|\s+/$1\|\t/' ids.txt
```

Then trim these down to complete genomes:

```
    egrep 'complete genome|complete sequence' ids.txt | grep -v plasmid | grep -v 'whole genome shotgun sequence' | grep -v NR_ > complete_genome_ids.txt
```

We edited this last file manually to ensure that we only have complete bacterial genomes.


We also downloaded the gbff files and orf files from RefSeq, and then used those to create a tbl file with both protein and CDS sequences, and create file called protein_list that has a list of all of the gbff files that we extracted.

```
    F=$(head -n $SGE_TASK_ID protein_list | tail -n 1)
    UF=$(echo $F | sed -e 's/.gz//')
    gunzip $F
    perl PhageHosts/code/genbank2flatfile.pl $UF
    gzip $UF
```

Based on this output we create two files, refseq_proteins.faa (with proteins) and refseq_orfs.faa (with DNA) from just the complete genomes.

```
    python PhageHosts/code/tbl2protdna.py .
```

## Extracting taxonomy information

We wrote code to automatically get the taxonomy for most of the hosts from the RefSeq files, but there were a few that we could not map, so we added those manually. To whit:

```
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
```

Then we generate the tsv file:

```
    python PhageHosts/code/phage_host_taxonomy.py  > phage_host_taxonomy.tsv.
```

We also just made two files with tuples of genome NC id and taxonomy id.

**NOTE:** *When we score the connection between phage and host, we need to know the taxonomy ID of the host, not that of the phage, to see if there is a match. Therefore, in this file the taxonomy ID is that of the **host** (not the phage!)*

```
    python PhageHosts/code/phage2taxonomy.py  > phage_host_taxid.txt
```

and another file refseq2taxonomy.py which added the tax id to the list of complete bacterial genomes, and then we made a list of bacteria and taxid:

```
    cut -f 2,4 /lustre/usr/data/NCBI/RefSeq/bacteria/complete_genome_ids_taxid.txt | perl -pe 's/^.*\|N/N/; s/\.\d+\|//' > bacteria_taxid.txt
```

We trimmed out any phages that we can not match at the species level:

```
    python2.7 PhageHosts/code/comparePhageToHosts.py
```

Then we combined those into a single file for all tax ids

```
    cat phage_taxid.txt bacteria_taxid.txt > data/all_host_taxid.txt
```

And add the taxonomy to those files:

```
    python PhageHosts/code/phage_host_add_taxonomy.py all_host_taxid.txt > data/all_host_taxid_taxonomy.txt
```
Note that in this process we deleted two phages whose hosts were not really known (NC_000935) APSE-1 whose host was Endosymbiont and Zamilon virophage (NC_022990) whose host is Mont1 megavirus)

### Resulting files

This resulted in one key file that we use in this work [data/all_host_taxid.txt](https://github.com/linsalrob/PhageHosts/blob/master/data/all_host_taxid.txt) which just has two columns, a RefSeq ID (which we sometimes call NC id) and a taxonomy ID. When the RefSeq ID refers to a bacterial sequence, the taxonomy ID is of the bacterium from where the sequence came. When the RefSeq ID refers to a phage sequence, the taxonomy ID is of the phage's host.

## Scoring predictions

In general, we are going to print out tab separated files containing the phage ID in the first column, and the ID of any potential hosts in subsequent columns. We do not necessarily know how many columns there will be. This is a flexible format that allows us to convert those IDs to taxonomy IDs, and then use the NCBI hierarchy to move through the phylogenetic tree to score how good our matches are.

We have a few keys pieces of code for this work:

* NC2taxid.py is python code that converts any RefSeq ID (NC_\d+) into its associated taxonomy ID based on the association above, the bacterial RefSeq IDs map to bacterial taxonomy IDs and the phage RefSeq IDs map to host taxonomy IDs.
* scoreTaxId.py is python code that takes the a set of taxonomy IDs and scores all subsequent elements in the set against the first member of the set, using 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom' taxonomic levels from NCBI.


# 1. Genetic Homology

## Comparing the phages against the complete bacterial genomes (blastx)

Start by making a database of just the complete genomes protein sequences

```
    PhageHosts/code/refseq2complete.py $HOME/phage/host_analysis/bacteria_taxid.txt  refseq_proteins.faa complete_genome_proteins.faa
```


and then blast those:

```
    PhageHosts/code/split_blast_queries_edwards_blastplus.pl -f ../../phage_with_host.fna -n 100 -p blastx -d phage.blastx -db RefSeq/bacteria/complete_genome_proteins.faa -evalue 0.01
    cat phage.blastx/*blastx > phage.complete.blastx
```

*NOTE:* There are mutliple proteins with the same ID that come from different genomes, however these are identical proteins (at the amino acid level). Therefore the blastx searches all return the same results for each protein. I use this one line of Perl code to print out unique solutions from the blastx output.

```
    perl -ne 'print unless ($s{$_}); $s{$_}=1' phage.complete.blastx > phage.complete.unique.blastx
```


We convert the output so we just have NC identifiers:

```
    python PhageHosts/code/blastx2genome.py phage.complete.unique.blastx phage.genomes.blastx
```


Now we just count the hosts with the most number of hits to each of the phages, and score those hits

```
    python PhageHosts/code/mostBlastHits.py phage.genomes.blastx > most.tax
    python2.7 PhageHosts/code/scoreTaxId.py most.tax > score.tax
```

## Comparing phages to the non-redundant (nr) database.

To see what would happen, we also compared the phages to all the proteins in the non-redundant database. Without going into the results in any detail, they weren't as good as using the complete genomes because there are more proteins in the nr database and thus more confusion from the similarity searches.

We used blast to compare the complete phage genomes against all bacterial proteins in the GenBank nr and then gi_taxid table of the [NCBI Taxonomy Site](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) ([here is an alternate link](http://tinyurl.com/ncbiftp) as there is an issue linking to FTP sites from github) to get the taxonomy id of the top hits. We compared those hits with the expected hosts.

The blast is a standard blast using nr as the database, and our output file is in blast tab separated output format (standard format with *qlen* and *slen* added to the output). We then use

```
    python PhageHosts/code/blast2taxid.py  phage_with_host.fna.blastx prot > phage.host.blastx.taxid
```

to create a list of phage and host taxids. Now we have to convert these to taxids and score the hits:

```
    for i in all equal best; 
	    do echo $i; 
    	python PhageHosts/code/blast_hits.py phage.host.blastx.taxid $i ${i}_hits.txt; 
	    python2.7 PhageHosts/code/scoreTaxId.py ${i}_hits.txt > ${i}_hits.score;
    done
```

This creates three files, `all_hits.score`, `equal_hits.score`, `best_hits.score`.


## Similarity to complete genomes (blastn)

We start with direct matches between the phages and hosts using the complete bacterial genomic DNA, first by making a database and then by using our cluster to blast the phages against the database. *Note:* We change the default values here to make them the similar as NCBI blast for a whole genomes, however that usually uses a word size of 28. I cut the word size to 5 so we get some shorter matches. Alternate parameters will affect the results!

```
    make a blast db: makeblastdb -in bacteria.genomic.fna -dbtype nucl
    PhageHosts/code/split_blast_queries_edwards_blastplus.pl -f phage_with_host.fna -n 100 -db /lustre/usr/data/NCBI/RefSeq/bacteria/complete_genomes.fna -d phage.blastn -p blastn -word_size 5 -evalue 10 -num_descriptions 100 -reward 1 -penalty -2 -gapopen 0 -gapextend 0 -outfmt '6 std qlen slen'
    cat phage.blastn/blastn >> phage_host.blastn
```

We then figure out the best hit for each phage, and score those hits.

```
    python2.7 PhageHosts/code/parse_blastn.py phage_host.blastn > phage.hits
    python PhageHosts/code/NC2taxid.py phage.hits > phage.taxid
    python2.7 PhageHosts/code/scoreTaxId.py phage.taxid > blastn.score
```

To make a table of all the scores between all phage genomes and all bacterial
genomes we use 

```
    python2.7 ../GitHubRepository/PhageHosts/code/blastn_all_v_all.py phage_host.blastn
```


# 2. CRISPR spacers. CRISPR spacers

Pilercr was downloaded from [drive5](http://www.drive5.com/pilercr/) and run against all complete genomes to create the database [data/pilercr1.06.output.fna](data/pilercr1.06.output.fna). 

blastn searches were performed using modified parameters that should be better for short matches:

```
    /usr/local/blast+/bin/blastn -db pilercr1.06.output.fna -query phage_with_host.fna -out phage_with_host.fna.crispr.blastn -outfmt '6 std qlen slen' -evalue 1 -gapopen 10 -penalty -1 -gapextend 2 -word_size 7 -dust no -task blastn-short
```

These blast searches were [compared](PhageHosts/code/score_blast.py score_blast.py), [converted to taxids](PhageHosts/code/crispr_blast2tax.py crispr_blast2tax.py), and [scored](PhageHosts/code/scoreTaxId.py scoreTaxId.py):

```
    python PhageHosts/code/score_blast.py phage_with_host.fna.crispr.blastn best > crispr.best.hits
    python PhageHosts/code/crispr_blast2tax.py bas.best.hits > crispr.taxid
    python2.7 PhageHosts/code/scoreTaxId.py bas.taxid  > bas.score
```


# 3. Exact matches

To identify the longest exact matches between the phage and bacteria we took a two step approach. This is easier to code than the complete solution, and although it takes slightly longer to run, the runtime is not an important concern in this analysis. 

We start by finding all 15mer hits between the bacteria and the phages, and compiling those in order of their location. We then combine those identical hits into longer stretches of similarity, because consecutive, overlapping, hits must have come from two longer matching sequences.

To find all 15-mers in common between the phage and bacteria, sort and join the matches, and list all the hits we use the following code. Note that we need to find the matches in both the forward and reverse direction, and then combine those matches.

```
    perl PhageHosts/code/find_exact_matches.pl ../phage_with_host.fna RefSeq/bacteria/complete_genomes.fna > phage.15mers.bacteria.txt
    perl PhageHosts/code/find_exact_matches.pl ../phage_with_host.rc.fna RefSeq/bacteria/complete_genomes.fna > phage.15mers.bacteria.rc.txt
    perl PhageHosts/code/sort_and_join_exact_matches.pl phage.15mers.bacteria.txt > phage.kmers.bacteria.txt
    perl PhageHosts/code/sort_and_join_exact_matches.pl phage.15mers.bacteria.rc.txt >> phage.kmers.bacteria.txt
    perl PhageHosts/code/longest_exact_match.pl phage.kmers.bacteria.txt > phage_best_hits.txt 
```

to find the longest match per phage, tabulate the RefSeq IDs, and score the matches we use:

```
    perl PhageHosts/code/longest_exact2tbl.pl  > phage_best_hits_NCIDS.txt
    python PhageHosts/code/NC2taxid.py phage_best_hits_NCIDS.txt  > phage_best_hits_taxa.txt
    python2.7 PhageHosts/code/scoreTaxId.py phage_best_hits_taxa.txt > score.txt
```

To create the supplementary figure, we used 

```
    python2.7 ../GitHubRepository/PhageHosts/code/exact_match_plot.py phage.kmers.bacteria.txt  > phage_kmers.counts
```

# 4. Oligonucleotide profiles

## GC Content of Coding Regions

Similar to codon usage below, we calculate the G+C/G+C+A+T content of the coding regions. This gives us a simple one dimensional matrix that we can use to map phages to hosts.

To count the GC content:

```
    python PhageHosts/code/cds_gc.py ../phage_with_host.cds.fna  > phage.gc.tsv
```

To count the bacterial GC content and limit that to complete bacterial sequences we used:

```
    python PhageHosts/code/cds_gc.py NCBI/RefSeq/bacteria/bacteria.$SGE_TASK_ID.orfs.fna.gz > bacteria/bacteria.$SGE_TASK_ID.gc.tsv
    python PhageHosts/code/complete_bacteria.py  > complete_bacteria_gc.tsv
```

Now we have two files: `complete_bacteria_gc.tsv` and `phage.gc.tsv` and we use pairwise [Euclidean distance](http://en.wikipedia.org/wiki/Euclidean_distance) to find the closest organisms, convert them to taxa, and score the results.

```
    python PhageHosts/code/cds_distance.py phage.gc.tsv complete_bacteria_gc.tsv  > closest_genomes.ids
    python PhageHosts/code/NC2taxid.py closest_genomes.ids > closest_genomes.taxa
    python PhageHosts/code/scoreTaxId.py closest_genomes.taxa > score.txt
```


## Codon usage

We count the codon profiles in both host and bacteria, and then calculate the [Euclidean distance](http://en.wikipedia.org/wiki/Euclidean_distance) using [codon_distance.py](code/codon_distance.py) to identify the closest host for each phage:

```
    python PhageHosts/code/codon_usage.py ../phage_with_host.cds.fna  > phage_codon_counts.tsv
    python codon_usage.py $DIR/bacteria.$SGE_TASK_ID.orfs.fna.gz > bacteria/bacteria.$SGE_TASK_ID.codons.tsv
    python PhageHosts/code/codon_distance.py phage_codon_counts.tsv bacteria_codon_counts.tsv  > phage_bacteria_predictions.out
```


convert that output to taxonomy IDs and score the output:

```
    python PhageHosts/code/NC2taxid.py phage_bacteria_predictions.out > phage_bacteria_predictions.taxid
    python2.7 PhageHosts/code/scoreTaxId.py phage_bacteria_predictions.taxid > score.txt
```


## Kmer profiles

We used [Jellyfish](http://www.genome.umd.edu/jellyfish.html) to count *k-*mers in the DNA sequences. The longest phage sequence is 358,663 bp, so we used a hash size of 400,000 to keep the whole thing (or most of it) in memory.

Initially we started with this script to list all the 15-mers:

```
    for f in $(ls ../phage_with_host.fna.files_NC/);
    	do n=$(echo $f | sed -e 's/\.\.\/phage_with_host.fna.files_NC\///; s/fna$/tsv/');
    	echo $f $n;
	    jellyfish count -s 400000 -t 32 -C -m 15 -o 15mers.txt $f;
	    jellyfish dump -ct 15mers.txt_0 > 15mers/$n;
    done
```

and then we modify it to count all k-mers between 3 and 20. That runs on a node of the cluster in a few minutes.


We then combine those outputs to a single file using [combineKmerCounts.py](code/combineKmerCounts.py):

```
    python2.7 PhageHosts/code/combineKmerCounts.py phage_kmer_counts/15mers 15mers.txt 0
```

That makes a single .txt file for each kmer size.

To count all bacterial kmers we need to extract the fasta sequences to individual files and run jellyfish on all of those. The combination of bash and perl does that:

```bash
    BACTGZ=$ENV{GENOME}
    BACT=$(echo $BACTGZ | sed -e 's/.gz//')

    if [ ! -e kmers ]; then mkdir kmers/; fi
    if [ ! -e sequences ]; then mkdir sequences/; fi

    mkdir sequences/
    ./split.pl $BACTGZ sequences/

    mkdir kmers/
    for FILE in $(ls sequences/); do
        for k in $(seq 3 10); do 
                jellyfish count -s 400000 -t 32 -C -m $k -o kmers/$FILE.${k}kmer sequences/$FILE
                if [ -e kmers/$FILE.${k}kmer_1 ]; then
                        jellyfish merge -o kmers/$FILE.${k}kmer.output.jf kmers/$FILE.${k}kmer_*
                else
                        mv kmers/$FILE.${k}kmer_0 kmers/$FILE.${k}kmer.output.jf
                fi
                jellyfish dump -ct kmers/$FILE.${k}kmer.output.jf > kmers/$FILE.${k}kmer.tsv
                rm -f kmers/$FILE.${k}kmer_*  kmers/$FILE.${k}kmer.output.jf 
        done
    done
```



This generates a single tsv file for every genome that we are going to use in the analysis.

Finally, we just make sure we only include genomes that we are interested in for our analysis.

```
    for i in $(ls all_phages); do echo $i; python trim.py all_phages/$i all_phages_trimmed/$i; done
    for i in $(seq 3 9); do python trim.py bacterial_genomes/kmers/$i.kmers bacterial_kmers/$i.kmers; done
```

# 5. Cooccurrence analysis


We need to predict the presence of bacteria and phage in all of the metagenomes that we have downloaded. We have 3,205 metagenomes that we have downloaded from [MG-RAST](http://metagenomics.anl.gov/), however we want to only use our complete set of phages and/or bacteria to do the predictions, and not the random set of organisms that are present in MG-RAST.

## Predicting bacterial presence

We start with [FOCUS](http://edwards.sdsu.edu/FOCUS/), but we need to make our unique databases first. Download the FOCUS code, and make the databases like this. The final FOCUS command makes the bacterial database.

```
    cp -r ~/bin/focus/bin focus_bacteria
    cd focus_bacteria/db/
    mv k6 k6_old
    mv k7 k7_old
    mv k8 k8_old
    head -n 1 k6_old > k6
    head -n 1 k7_old > k7
    head -n 1 k8_old > k8
    cd ../..
    grep -r \> NCBI/RefSeq/bacteria/complete_genomes.fna.files.NC/  > complete_genomes_files.txt
    perl PhageHosts/code/focus_grep2list.pl < complete_genomes_files.txt  > complete_genomes_files.focus.txt
    python2.7 focus_bacteria/focus.py -d complete_genomes_files.focus.txt
```

Now we can run all the bacterial genomes using focus. This script makes a directory for each job, cd's there, copies the data, runs focus, and then cleans up.

This is focus.sh

```bash
    #!/bin/bash

    export LD_LIBRARY_PATH=:$HOME/lib
    GZFILE=$(head -n  $SGE_TASK_ID mg-rast.txt | tail -n 1)
    mkdir /scratch/$SGE_TASK_ID
    cd /scratch/$SGE_TASK_ID
    gunzip -c $GZFILE > $SGE_TASK_ID.fasta
    echo "Input: " $GZFILE  > cooccurence/focus_NC/$SGE_TASK_ID.focus
    python2.7 focus.py -m 0 -q $SGE_TASK_ID.fasta >> cooccurence/focus_NC/$SGE_TASK_ID.focus
    cd $HOME/phage/host_analysis/cooccurence
    rm -rf /scratch/$SGE_TASK_ID
```

Since there are 3025 files in mg-rast.txt, we use qsub to run this command `qsub  -cwd -S /bin/bash -t 1-3025:1 -o sge -e sge ./focus.sh`

This takes just a few hours to predict all the bacteria present in all the metagenomes that have been produced, and once the computation is complete, we combine them into a single file with the bacterial predictions normalized to the lengths of the metagenomes: `python ~/bioinformatics/phage_host/focus_parse.py focus_NC2 > focus_bacteria_predictions.tsv`

## Predicting phage sequences

We have already run megablast for all phages against all metagenomes for another project, and so we just use this as our input.

The code below makes two tables called `raw_phage_mg_counts.tsv` and `normalized_phage_mg_counts.tsv` listing the presence of phages in each of the metagenomes. The first is just the number of times that phage is seen in each metagenome, and the second is the number divided by the number of reads in the metagenome. It is the normalized number we use (longer metagenomes have more hits, after all).

```
    python2.7 PhageHosts/code/count_phage_mg_hits.py 2> err
```


Our computations produce two files, `normalized_phage_mg_counts.tsv` and `focus_bacteria_predictions.tsv` that we need to compare. To do so, we calculate Pearson's correlation between genomes and phages, convert those to taxonomy IDs and then score them.

```
    python PhageHosts/code/scoreDistances.py pearson_correlations.tsv max > pearson.ncids
    python PhageHosts/code/NC2taxid.py pearson.ncids > pearson.tids
    python2.7 PhageHosts/code/scoreTaxId.py pearson.tids > pearson.score
```

We have also developed code to compare the scores using Euclidean distance, Bray Curtis distance, City Block distance, Hamming distance, and Jaccard distance, and we can run those sequentially, performing the calculations, computing the taxonomy IDs and scoring the results.

```
    for i in euclidean braycurtis cityblock hamming jaccard; 
	    do python PhageHosts/code/scoreDistances.py ${i}_correlations.tsv min > ${i}_correlations.ncids;
    	python PhageHosts/code/NC2taxid.py ${i}_correlations.ncids > ${i}_correlations.tids;
    	python2.7 PhageHosts/code/scoreTaxId.py ${i}_correlations.tids > ${i}_correlations.score;
    done
```


## Mutual Information

This is based on discussions with Peter Salamon, who said:

"You need a threshold of present in the metagenome to turn your current variables into 0,1 variables that signal presence of the entity in the metagenome. As a start, you could threshold on mean value for that entity (phage or bacterium). Then you need the fraction of metagenomes in which the phage and the bacterium co-occur, the fraction with one not the other and the fraction with the other and not the one and the fraction with neither. This gives the joint distribution of the four cases. You also need the marginal distributions of fraction in which the phage occurs and the fraction in which the bacterium occurs. Then simply calculate the Kullback entropy of the joint distribution relative to the product distribution."

```
    python2.7 PhageHosts/code/mutual_information.py normalized_phage_mg_counts_scaled.tsv focus_bacteria_predictions.tsv 0 > mi.mean.ncids
    python PhageHosts/code/NC2taxid.py mi.mean.ncids > mi.mean.tids
    python2.7 PhageHosts/code/scoreTaxId.py mi.mean.tids > mi.mean.score
```
