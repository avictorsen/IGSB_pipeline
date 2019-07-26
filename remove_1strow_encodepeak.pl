#!/usr/bin/perl
BEGIN {
	push @INC, "/home/malijia/bin";
}
use strict;
use Lijia;

# ljma@uchicago.edu
# 
#  

if (@ARGV != 2) {
	print "usage: .pl InDir OutDir\n";
	exit;
}

my $indir = shift;
my @Peaks=`find $indir -maxdepth 1 -name "*.encodePeak"`;
my $outdir = shift;

Ptime("start!");

foreach my $file (@Peaks) {
	chomp($file);
#	print "remove first line of $file\n";
	open(IN,$file) or die $!;
	open(OUT,">$outdir/$file.tmp") or die $!;
	while (<IN>) {
		if($_=~/^track/){
			my $removed=$_;
			print "$file\n";
			print "removed:\t$removed";
		}else{
			print OUT "$_";
		}
	}
	close IN;
	`mv $file.tmp $file\n`;
}

Ptime("Done!");


