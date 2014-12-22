#!/usr/bin/perl -w

# get everything out of a genbank file

use Bio::SeqIO;
use strict;
use Rob;
my $rob = new Rob;

my $usage=<<EOF;
$0 <list of genbankfiles>

EOF

die $usage unless ($ARGV[0]);


# get a list of sequences we need
my $cf='/usr/lustrefs/anthill/usr/data/NCBI/RefSeq/bacteria/complete_genome_ids.txt';
open(IN, $cf) || die "Can't open $cf";
my %want;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	map {$want{$_}=1} split /\|/, $a[1];
}
close IN;

my $c;
foreach my $file (@ARGV)
{
	my $outF = $file;
	$outF =~ s/gbff/tbl/;
	# read the cds from this sequence file
	my $fafile = $file;
	$fafile =~ s/genomic.gbff/orfs.fna.gz/;
	my $cds = $rob->read_fasta($fafile);
	map {my $s=$_; $s=~s/\s+.*$//; $cds->{$s}=$cds->{$_}} keys %$cds;
	if (-e $outF) {die "$outF already exists. Not overwriting\n"}
	open(TB, ">$outF") || die "can't open tbl";
	my $sio=Bio::SeqIO->new(-file=>$file, -format=>'genbank');
	while (my $seq=$sio->next_seq) {
		my $seqname=$seq->display_name;
		# $seqname is the refseq id and primary_id is the gi
		next unless ($want{$seq->primary_id()} || $want{$seqname});

		print STDERR "Parsing $seqname with id ",  $seq->primary_id(), "\n";
		my @seqdata = ( $seq->display_name.".".$seq->seq_version(), $seq->length(), $seq->desc(), $seq->primary_id());
		my $getDna=1;
		if (length($seq->seq) < 50000) {$getDna=0; print STDERR "Not able to get the DNA sequence from $seqname in $file\n"}
		foreach my $feature ($seq->top_SeqFeatures()) {
			
			$c++;
			my $id; # what we will call the sequence
			my ($trans, $gi, $geneid, $prod, $np);
			my $locus = "";

			eval {$trans = join " ", $feature->each_tag_value("translation")};
			eval {$np = join " ", $feature->each_tag_value("protein_id")};
			if (!$np) {
				eval {$np = join " ", $feature->each_tag_value("product")};
			}
			if ($trans && !$np) {print STDERR "No NP for $trans. Skipped\n"; next}
			elsif (!$trans && $np) {print STDERR "No translation for $np. Skipped\n"; next}
			next unless ($trans && $np);

			eval {
				foreach my $xr ($feature->each_tag_value("db_xref")) 
				{
					($xr =~ /GI/) ? ($gi = $xr) : 1;
					($xr =~ /GeneID/) ? ($geneid = $xr) : 1;
				}
			};

			eval {$locus = join " ", $feature->each_tag_value("locus_tag")};
			eval {$prod  = join " ", $feature->each_tag_value("product")};

			my $oids = "locus:$locus;"; 
			($geneid)  && ($oids.="$geneid;");
			($gi)      && ($oids.="$gi;");
			$oids =~ s/\;$//;

## columns: sequence name, version, length of sequence, genome name, genome id, protein id, protein start, protein stop, protein strand, protein sequence, protein function, protein aliases

			unless ($prod)  {print STDERR "No product for $np\n"; $prod="hypothetical protein"}
			#print STDERR "Getting sequence from ", $feature->start, " to ", $feature->end, " from seq length ", $seq->length(), "\n";
			my $dna="";
			if ($cds->{$locus}) {$dna=$cds->{$locus}}
			elsif ($getDna) {
				$dna = $seq->subseq($feature->start, $feature->end);
				if ($feature->strand == -1) {$dna=$rob->rc($dna)}
			}
			else {
				print STDERR "Could not find DNA for $seqname (locus:$locus)\n";
			}
			print TB join("\t", $seqname, @seqdata, $np, $feature->start, $feature->end, $feature->strand, $trans, $dna, $prod, $oids), "\n";
			# print TB join("\t", $feature->start, $feature->end, $feature->strand, $dna), "\n\n";
		}
	}
}


