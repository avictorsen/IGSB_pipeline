#!/usr/bin/perl
use strict;
#use SVG;

# ljma@uchicago.edu
# Oct-14-2010
# SAM statistics


if (@ARGV !=2) {
	print ".pl ./(BWA.SAM dir) out \n";
	exit;
}

my $infile = shift;
my @Files = `find $infile -name "*sam"`;
#open(MAP, "$infile") or die $!;
my $infile = shift;
open(OUT, ">>$infile") or die $!;

my $time = localtime;
print "started generating SAM statistics @ $time\n";

for my $file (sort {$a cmp $b} @Files) {
	chomp($file);
	open(MAP, "$file") or die $!;
	print OUT "$file\n";
	print "running $file\n";

	my($un,$multi,$multi_d,$unique,$unique_d,$unique_e,$unique_1,$unique_2,$unique_3)=(0,0,0,0,0,0,0,0,0);
	my(%Distinct_M,%Distinct_U);
	while (<MAP>) {
		chomp;
		my $line = $_;
		my @info = split/\t/,$line;
		next if ($info[0] =~ /^\@/);
		if ($info[2] eq "*") {
			$un ++;
		}elsif ($line=~/\tXA\:/) {
			$multi ++;
			if (($info[1]&16)==16) {
#			print "$info[0]\t$info[1]\n";
#			print "$info[9]\n";
				$info[9]=~tr/ATGCatgc/TACGtacg/;
				$info[9]=reverse($info[9]);
#			print "$info[9]\n";
			}else{}
			if (!exists $Distinct_M{$info[9]}) {
				$multi_d ++;
				$Distinct_M{$info[9]}=1;
			}else{}
		}else {
			$unique ++;
			my($strand,$mismatch_tag);
			if (($info[1]&16)==0) {
				$strand="-";
			}else{
				$strand="+";
			}
			if ($line=~/\t(MD\:Z\:\w+)/) {
				$mismatch_tag=$1;
			}
			my $alledit_tag="$strand-$info[2]-$info[3]-$info[5]-$mismatch_tag";  ##### strand-chr-start-indel-mismatch
				if (!exists $Distinct_U{$alledit_tag}) {
					$unique_d ++;
					$Distinct_U{$alledit_tag} = 1;
				}else{}

			if ($line=~/\tXM\:i\:(\d+)/) {
				if ($1 == 0) {
					$unique_e ++;
				}elsif ($1 == 1) {
					$unique_1 ++;
				}elsif ($1 == 2) {
					$unique_2 ++;
				}elsif ($1 == 3) {
					$unique_3 ++;
				}else{}
			}
		}
	}
	close MAP;

	my$total=$un+$multi+$unique;
	my$total_d=$multi_d+$unique_d;
	print OUT "Total\t$total\nmapable_distinct\t$total_d\nunmappable\t$un\nunique_map\t$unique\nunique_distinct\t$unique_d\nunique_exact\t$unique_e\nunique_1\t$unique_1\nunique_2\t$unique_2\nunique_3\t$unique_3\nmultiple_map\t$multi\nmultiple_distinct\t$multi_d\n";
        `rm -f $file`;
}

my $time = localtime;
print "finished @ $time\n";

close OUT;
