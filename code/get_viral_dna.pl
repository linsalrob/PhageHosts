#!/usr/bin/perl -w
#

use strict;
use Rob;
my $rob=new Rob;
open(IN, "refseq.txt") || die $!;
my %want;
my %phage;
my %euk;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	$want{$a[0]}=$a[2]  if ($a[2]);
	if ($a[4] eq "YES") {$phage{$a[0]}=1}
	elsif ($a[4] eq "NO") {$euk{$a[0]}=1}
	else {die "What is $a[4]!"}
}
close IN;

my %found;
my $fa=$rob->read_fasta("2014_07_21/refseq/viral.1.1.genomic.fna");
open(PHG, ">phage_with_host.fna") || die $!;
open(EUK, ">eukaryotic_with_host.fna") || die $!;

foreach my $i (keys %$fa) {
	my $k=0;
	map {s/\.\d+$//; $k=$_ if ($want{$_})} split /\|/, $i;
	next unless ($k);
	$found{$k}=1;
	if ($phage{$k}) {
	print PHG ">$i\n", $fa->{$i}, "\n" if ($k);
} elsif ($euk{$k}) {
	print EUK ">$i\n", $fa->{$i}, "\n" if ($k);
}
else {
	print STDERR "HUH?? $k\n";
}
	
}

foreach my $k (keys %want) {
	print STDERR "NOT FOUND: $k\n" unless ($found{$k});
}

