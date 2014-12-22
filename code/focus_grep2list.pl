while (<>) {
	chomp;
	s/^(.*?)\:/'$1'\t/;
	s/\>gi.*?(NC_\d+).*/$1/;
	@a=split /\t/;
	print $a[0] . "\t" . (($a[1] . "\t") x 8) . "\n";
}

