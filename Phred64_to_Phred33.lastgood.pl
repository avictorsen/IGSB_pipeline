#!/usr/bin/perl -w
use strict;
my $count = 1;
while(<>) {
	my $line = $_;
	if($count % 4 == 0) {# only look at every 4th line
	if ($line =~ /[J-h]/){
		$line =~ s/[(\;\<\>\=\?]/\!/g; # incase solexa
		$line =~ tr/@-h/!-I/; # for converting to Phred+33
		}
	}
	print $line;
	$count++;
}

