#!/usr/bin/perl 
#

use strict;
use Rob;

my $q=shift || die "$0 <query sequence> <database sequences> USE THE PHAGE FOR THIS (shorter genomes)";
my $d=shift || die "$0 <query sequence> <database sequences>";

#my $kmer = 15;
my $kmer = 11;

my $rob=new Rob;

my $fa = $rob->read_fasta($q);
{
	my $temp;
	map {my $i = $_; $i =~ s/^.*(NC_\d+).*$/$1/; $temp->{$i}=$fa->{$_}} keys %$fa;
	$fa = $temp;
}

my $oligo;
foreach my $id (keys %$fa) {
	my $co = $rob->oligos_in($fa->{$id}, $kmer);
	map {if (index($_, "N") == -1) {push @{$oligo->{$_}}, $id}} keys %$co;
}
$fa={};

{
	my @os = keys %$oligo;
	#print STDERR "OLIGOS\n=======\n", join("\n", @os[0 .. 10]), "\n";
}

## we're going to stream the fasta file and print out all matches for a sequence
open(IN, $d) || die "Can't open $d";
my $seq="";
my $id;
while (<IN>) {
	chomp;
	if (s/^>//) {
		if ($seq) {&printmatches($id, $seq)}
		$seq="";
		$id=$_;
		$id =~ s/^.*(NC_\d+).*$/$1/;
	} else {
		$seq .= uc($_);
	}
}
&printmatches($id, $seq);
close IN;

sub printmatches {
	my ($id, $seq)=@_;
	print STDERR "Scanning $id\n";
	my $posn=0;
	while ($posn + $kmer < length($seq)) {
		my $ss=substr($seq, $posn, $kmer);
		next if (index($ss, "N") > -1);
		if ($oligo->{$ss}) {
			map {print join("\t", $_, $id, $posn, $ss), "\n"} @{$oligo->{$ss}}
		}
		$posn++;
	}
}





