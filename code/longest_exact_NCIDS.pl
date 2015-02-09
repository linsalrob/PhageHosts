use strict;
my %w;
open(IN, "/home3/redwards/phage/host_analysis/phage_with_host.tsv") || die $!;
while (<IN>) {
	my @F=split /\t/;
	$w{$F[0]}=1;
}
close IN;

open(IN, "phage_best_hits.txt") || die $!;
while (<IN>) {
	my @F=split /\t/;
	$F[0] =~ s/.*(NC_\d+).*/$1/;
	$F[1] =~ s/.*(NC_\d+).*/$1/;
	print join("\t", @F);
}
close IN;

