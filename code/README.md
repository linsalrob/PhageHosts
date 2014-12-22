Code Used In the Analysis of Phages and their Hosts
===================================================

Almost all of this code was written by Rob Edwards or Bas Dutilh. The good stuff was written by Bas, the crappy stuff by someone else.

You are welcome to use this code in any way you wish, but we are not liable if it screws up your system. I don't think we have an rm -rf anywhere, but you never know!


0. General utility code
=======================

* check_dna.py -- check the DNA that we get from RefSeq
* comparePhageToHosts.py -- compare phage and their hosts (not used in the end)
* complete_bacteria.py -- identify all complete bacteria in data downloaded from RefSeq
* genbank2flatfile.pl
* missing_phage_blastn.py
* phage.py
* phage_with_specified_hosts.py
* phage2taxonomy.py
* phage_host_add_taxonomy.py
* phage_host_taxonomy.py
* plot.py
* refseq2complete.py
* refseq2taxonomy.py
* score_blast.py
* split_orfs_by_genome.py


* split_blast_queries_edwards_blastplus.pl -- split a DNA sequence file and run blast+ on the cluster
* tbl2protdna.py -- convert NCBI tbl format to proteins and DNA
* translate.pl -- translate an ORF

0a. Scoring code
----------------
This code is used to convert RefSeq NC_ ids to taxonomy ids and then to score the code
* NC2genus_species.py -- find the genus and species for our RefSeq IDs
* NC2taxid.py -- convert a list of RefSeq IDs to taxonomy IDs
* scoreDistances.py -- score the distances between genomes
* score.py -- legacy code - don't use this
* scoreTaxId.py -- convert the taxonomy IDs into a score
* score_taxonomy.py -- legacy code - don't use this


1. Protein similarity
=====================

1a. Similarity to known proteins (blastx)
-----------------------------------------

* blast2taxid.py -- convert blast results to taxonomy IDs
* blast2tax.py -- convert blast results to taxonomy IDs
* blast_hits.py -- get the hits out of blast files
* blastx2genome.py -- convert blastx results (eg against NR) to genomes
* mostBlastHits.py


1b.  Similarity to complete genomes (blastp)
--------------------------------------------

See the list above, since almost all the blast code is identical

2. Nucleotide Similarity (blastn)
==========================================

* parse_blastn.py
* tabulate_blastn.py
Also see the list above, since almost all the blast code is identical


3. Exact matches
======================================

* find_exact_matches.pl
* longest_exact2tbl.pl
* longest_exact_match.pl
* sort_and_join_exact_matches.pl

4. CRISPR Sequences
===================

* crispr_blast2tax.py
* score_crispr.pl

5. GC Content of Coding Region
===============================

* cds_gc.py -- calculate the cds %GC

6. Codon usage
==============
* cds_distance.py -- calculate the distance between cds codon usage
* codon_distance_one.py -- calculate a single codon distance
* codon_distance.py -- calculate the distance between all codon usages
* codon_usage.py -- calcualte the codon usage

7. Kmer profiles
================

* combineKmerCountsBacteria.py -- combine all the kmer counts from JellyFish into a single table
* combineKmerCounts.py -- combine all the kmer counts from JellyFish into a single table

8. Cooccurrence analysis
=========================

* cooccurrence_all.py
* cooccurrence_pearson.py
* count_phage_mg_hits.py
* extract_phage_mg_blast.py
* focus_grep2list.pl
* focus_parse.py
* focus_summary.py
* mutual_information.py


11. Methods that we did not use
===============================

These programs were written for some additional methods, but we did not use them in the paper.

11a. Codon adaptation index
---------------------------
* best_cai.py -- codon adaptation index. This was not used in the end
* cai.py -- calculate the codon adaptation index
* choose_genomes_with_highly_expressed_proteins.py -- find some highly expressed proteins for the CAI

11b. frac_q and frac_d
----------------------
* fracDQ_best.py
* fracDQ.py



