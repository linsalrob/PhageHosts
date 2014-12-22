#!/usr/bin/perl 
#
# score the crispr blast hits. We will use fractional counts for this
# for example: for i in $(grep NC_019422.1 crispr.blastn.blastn | cut -f 2); do grep -w $i /home/db/CRISPR/crispr.u-psud.fr/id.map.gs ; done
# has this output:
#	 5949	NC_017174_6_17	
#	 60301	NC_017174_12_16	
#	 41839	NC_013315_15_10|NC_013316_15_10|NC_017178_15_10|NC_017179_15_10	Clostridium difficile
#	 3662	NC_013974_10_5	Clostridium difficile
#	 46173	NC_017174_12_19	
#	 21655	NC_017174_12_21	
#	 7152	NC_013974_13_3	Clostridium difficile
#	 62857	NC_017173_10_2	
#	 52166	NC_017177_13_5	Clostridium difficile
#	 56039	NC_017174_6_31	
#
# We will ignore any blanks (these are either Archaea or more likely deleted entries from RefSeq), and then sum up what is left
use strict;
use Getopt::Std;
my %opts = ('p' => 0, 'f'=>0);
getopts('b:p:f:s', \%opts);

unless ($opts{b}) {
	die <<EOF;
$0 
-b crispr blastn output file
-p minimum percent identity of the hit (>=100)
-f minimum fractional coverage of alignment length/subject length (>=1)
-s just print out the summary
EOF
}

if ($opts{f} > 1) {
	print STDERR "Corrected -f $opts{f} to ", $opts{f}/100, "\n";
	$opts{f} /= 100;
}

my %host;
open(IN, "phage_with_host.tsv") || die $!;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	$host{$a[0]}=$a[2];
}
close IN; 

my %crispr;
open(IN, "/home/db/CRISPR/crispr.u-psud.fr/id.map.gs") || die $!;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	next unless ($a[2]);
	my @hosts = split /\|/,$a[2];
	$crispr{$a[0]}=\@hosts;
}
close IN; 

my %hostsfor; my %hits;
open(IN, $opts{b}) || die $!;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	next if ($a[2] < $opts{p});
	unless ($#a == 13) {die "You must include qlen and slen in the blast output"}

	next if (($a[3]/$a[13]) < $opts{f});
	$a[0] =~ /(NC_\d+)/;
	my $id = $1;
	map {if ($_) {$hostsfor{$id}{$_}++; $hits{$id}++}} @{$crispr{$a[1]}}
}	
close IN;

my ($same, $diff)=(0,0);
foreach my $id (keys %hostsfor) {
	$opts{s} || print "$id\t$host{$id}";
	my $best;
	foreach my $h (sort {$hostsfor{$id}{$b} <=> $hostsfor{$id}{$a}} keys %{$hostsfor{$id}}) {
		if ($best && $hostsfor{$id}{$best} == $hostsfor{$id}{$h} && $h eq $host{$id}) {$best=$h}
		unless ($best) {$best = $h}
		$opts{s} || print "\t$h (", (($hostsfor{$id}{$h}/$hits{$id})*100), " %)";
	}
	$opts{s} || print "\n";
	if ($best eq $host{$id}) {$same++}
	else {$diff++}
}

print STDERR "Correct assertions: $same Wrong assertions: $diff Total assertions: ", ($same+$diff), "\n";
