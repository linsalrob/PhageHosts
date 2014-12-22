#!/usr/bin/perl -w
#
# we need to sort the exact oligo matches and join up neighboring matches to make them longer
#
# The matches are ordered:
# [phage, genome, position in genome, sequence]
# we will read in all the hits for a genome, and then process those (to save memory). We sort them first by phage then by position

use strict;
use Data::Dumper;

my $f=shift || die "$0 <file of kmer hits to merge>";
open(IN, $f) || die $!;
my $lastgenome="";
my $hit;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	if ($lastgenome eq $a[1]) {
		push @{$hit}, \@a;
	}
	else {
		&processhits($hit) if ($hit);
		$lastgenome=$a[1];
		$hit = [\@a];
	}
}
&processhits($hit);

sub processhits {
	my $h = shift;
	my @sorted = sort {$a->[0] cmp $b->[0] || $a->[2] <=> $b->[2]} @$h;
	#print Dumper(\@sorted);
	my $i=0;
	my %skip;
	while ($i < $#sorted) {
		my $p=$sorted[$i]->[2];
		my $s=$i;
		while ($sorted[$i+1]->[2] == $p+1) {
			$p++;
			$i++;
			$skip{$i}=1;
			$sorted[$s]->[3] .= substr($sorted[$i]->[3], -1);
		}
		push @{$sorted[$s]}, length($sorted[$s]->[3]);
		print join("\t", @{$sorted[$s]}), "\n";
		$i++;
	}
}




