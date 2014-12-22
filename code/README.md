Code Used In the Analysis of Phages and their Hosts
===================================================

Almost all of this code was written by Rob Edwards or Bas Dutilh. The good stuff was written by Bas, the crappy stuff by someone else.

You are welcome to use this code in any way you wish, but we are not liable if it screws up your system. I don't think we have an rm -rf anywhere, but you never know!


0. General utility code
=======================

* check_dna.py -- check the DNA that we get from RefSeq
* comparePhageToHosts.py -- compare phage and their hosts (not used in the end)
* complete_bacteria.py -- identify all complete bacteria in data downloaded from RefSeq

* split_blast_queries_edwards_blastplus.pl -- split a DNA sequence file and run blast+ on the cluster
* tbl2protdna.py -- convert NCBI tbl format to proteins and DNA
* translate.pl -- translate an ORF

1. Protein similarity
=====================

1a. Similarity to known proteins (blastx)
-----------------------------------------

* blast2taxid.py -- convert blast results to taxonomy IDs
* blast2tax.py -- convert blast results to taxonomy IDs
* blast_hits.py -- get the hits out of blast files
* blastx2genome.py -- convert blastx results (eg against NR) to genomes

1b.  Similarity to complete genomes (blastp)
--------------------------------------------



2. Nucleotide Similarity (blastn)
==========================================


3. Exact matches
======================================


4. CRISPR Sequences
===================




5. GC Content of Coding Region
===============================

* cds_gc.py -- calculate the cds %GC



6. Codon usage
==============
* cds_distance.py -- calculate the distance between cds codon usage
* codon_distance_one.py -- calculate a single codon distance
* codon_distance.py -- calculate the distance between all codon usages
* codon_usage.py -- calcualte the codon usage





7. Codon adaption index (CAI)
=============================



8. Kmer profiles
================

* combineKmerCountsBacteria.py -- combine all the kmer counts from JellyFish into a single table
* combineKmerCounts.py -- combine all the kmer counts from JellyFish into a single table



9. Frac_D and Frac_Q
===================



10. Cooccurrence analysis
=========================




11. Methods that we did not use
===============================

These programs were written for some additional methods, but we did not use them in the paper.

11a. Codon adaptation index
---------------------------
* best_cai.py -- codon adaptation index. This was not used in the end
* cai.py -- calculate the codon adaptation index
* choose_genomes_with_highly_expressed_proteins.py -- find some highly expressed proteins for the CAI




Codon

* cooccurrence_all.py
* cooccurrence_pearson.py
* count_phage_mg_hits.py
* crispr_blast2tax.py
* extract_phage_mg_blast.py
* find_exact_matches.pl
* focus_grep2list.pl
* focus_parse.py
* focus_summary.py
* fracDQ_best.py
* fracDQ.py
* genbank2flatfile.pl
* longest_exact2tbl.pl
* longest_exact_match.pl
* missing_phage_blastn.py
* mostBlastHits.py
* mutual_information.py
* NC2genus_species.py
* NC2taxid.py
* parse_blastn.py
* phage2taxonomy.py
* phage_host_add_taxonomy.py
* phage_host_taxonomy.py
* phage.py
* phage_with_specified_hosts.py
* plot.py
* README.md
* refseq2complete.py
* refseq2taxonomy.py
* score_blast.py
* score_crispr.pl
* scoreDistances.py
* score.py
* scoreTaxId.py
* score_taxonomy.py
* sort_and_join_exact_matches.pl
* split_orfs_by_genome.py
* tabulate_blastn.py

