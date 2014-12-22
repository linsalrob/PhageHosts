use lib '/home3/redwards/bioinformatics/Modules';
use raeseqlib qw/translate_seq/;
use Rob;
my $rob = new Rob;
my $rae = new raeseqlib;
my $faf = shift || die "fasta file?";
my $fa=$rob->read_fasta($faf);
foreach my $i (keys %$fa) {
	my $t=translate_seq($fa->{$i}, 1,);
	$t =~ s/\*$//;
	print ">$i\n", $t, "\n";
}

