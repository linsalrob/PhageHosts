#!/usr/bin/perl -w
#
# TAke the output from ~/bioinformatics/phage_host/sort_and_join_exact_matches.pl and just list the longest match for each phage. If there are multiple equal longest matches, keep them all
#

use strict;
my $f = shift || die "$0 <output from sort_and_join_exact_matches.pl>";
open(IN, $f) || die $!;
my %len;
my %long;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	unless (defined $len{$a[0]}) {$len{$a[0]}=0}
	if ($a[$#a] > $len{$a[0]}) {
		$len{$a[0]}=$a[$#a];
		@{$long{$a[0]}} = ([$a[1], $a[$#a]]);
	}
	elsif ($a[$#a] == $len{$a[0]}) {
		push @{$long{$a[0]}}, [$a[1], $a[$#a]];
	}
}
foreach my $phage (keys %len) {
	foreach my $res (@{$long{$phage}}) {
		print join("\t", $phage, @$res), "\n";
	}
}


