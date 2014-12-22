#!/usr/bin/perl
#
use strict;


open(IN, "phage_best_hits.txt") || die $_;
my %h;
while(<IN>) {
	chomp;
	my @a=split /\t/;
	$a[0] =~ s/^.*ref\|([A-Z]+\_\d+).*$/$1/;
	$a[1] =~ s/^.*ref\|([A-Z]+\_\d+).*$/$1/;
	$h{$a[0]}{$a[1]}=1;
}

foreach my $k (keys %h) {
	print join("\t", $k, keys %{$h{$k}}), "\n";
}
