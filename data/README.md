Data Sets for the Analysis of Phages and their Hosts
====================================================

**NOTE:** In this work we are trying to find the right host for a phage. The way that we do that is predict the phage host and then for each phage we generate a list of RefSeq ID(s) of potential hosts. We convert that list to taxonomy IDs. When we do this we list the taxonomy ID of the host, not the taxonomy ID of the phage, and that way we can compare our predictions.

The following data were used in this analysis
* all_host_taxid_taxonomy.txt -- a list of all the bacteria, and the phages, and the taxonomy strings. The phage host in the second column is the taxonomy id of the host of the phage in the first colum. This file contains both bacteria and phages mixed together.
* all_host_taxid.txt -- a list of all the bacteria, and the phages, and the taxonomy strings. The phage host in the second column is the taxonomy id of the host of the phage in the first colum. This file contains both bacteria and phages mixed together.
* bacteria_taxid.txt -- a list of the bacteria that we use the in the analysis and the taxonomy id of those bacteria
* phage_host_taxid_allpresent.txt -- a list of the phages that have a genome in the data set and the taxonomy id of the host for the phage
* phage_with_host.tsv -- a list of the phages and the hosts that they infect
