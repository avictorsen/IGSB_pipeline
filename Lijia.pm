package Lijia;

use strict;

require Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 1.00;
@ISA = qw(Exporter);
@EXPORT = qw(Ptime Overlap Merge Quantile median avg sum min max roundup);
@EXPORT_OK = qw(Ptime Overlap Merge Quantile median avg sum min max roundup);


sub Ptime {
	my $time = localtime;
	my ($msg) = @_;
	print "$msg at $time\n";
}


sub Overlap{
#########out format "b which b-a",  "b wich aIb";  Intersection
	my ($as,$ae,$bs,$be) = @_;
	my %Data;
	for(my $i=0;$i<@$as;$i++){
		if($$as[$i]>$$ae[$i]){
			print "Error,a:start>end:$$as[$i]>$$ae[$i]\n";
			next;
		}
		for(my$j=$$as[$i];$j<=$$ae[$i];$j++){
			$Data{$j} .= "$$as[$i]-$$ae[$i];";
		}
	}
	my %Overlap;
	for(my $i=0;$i<@$bs;$i++){
		if($$bs[$i]>$$be[$i]){
			print "Error,a:start>end:$$bs[$i]>$$be[$i]\n";
			next;
		}
		for(my$j=$$bs[$i];$j<=$$be[$i];$j++){
			my @peaks = split/;/,$Data{$j};
			foreach my$p (@peaks){
				$Overlap{"$$bs[$i]-$$be[$i]"}{$p} = 1;
			}
		}
	}
	my %b_a;
	my %bIa;
	for(my $i=0;$i<@$bs;$i++){
		my @peaks =sort{my @a=split/-/,$a;
				my @b=split/-/,$b;
				if($a[0] == $b[0]){
					$a[1] <=> $b[1];
				}else{
					$a[0] <=> $b[0];
				}
			} keys %{$Overlap{"$$bs[$i]-$$be[$i]"}};
		if(@peaks == 0){
			$b_a{"$$bs[$i]-$$be[$i]"} = 1;
		}else{
			$bIa{"$$bs[$i]-$$be[$i]"} = \@peaks;
		}
	}
	return (\%b_a,\%bIa);
}

sub Merge{
#########out format "a-b", "b-a", "aUb", "aIb"; Union, Intersection
	my ($as,$ae,$bs,$be,$out) = @_;
	my %Data;
	my $j=0;
	for(my $i=0;$i<@$as;$i++){
		if($$as[$i]>$$ae[$i]){
			print "Error,a:start>end:$$as[$i]>$$ae[$i]\n";
			next;
		}
		$Data{$$as[$i]} .= "a_$j;";
		$Data{$$ae[$i]} .= "b_$j;";
		$j++;
	}
	for(my $i=0;$i<@$bs;$i++){
		if($$bs[$i]>$$be[$i]){
			print "Error,a:start>end:$$bs[$i]>$$be[$i]\n";
			next;
		}
		$Data{$$bs[$i]} .= "c_$j;";
		$Data{$$be[$i]} .= "d_$j;";
		$j++;
	}
	my $flaga=0;
	my $flagb=0;
	my %DataS;
	foreach my $d (sort{$a<=>$b} keys %Data){
#               print "raw:$d\t$Data{$d}\n";
		my @info = split/_|;/,$Data{$d};
		my $tempa=0;
		my $tempb=0;
		for(my$i=0;$i<@info;$i+=2){
			if($info[$i] eq "a"){
				$flaga++;
				$tempa=1;
			}elsif($info[$i] eq "b"){
				$flaga--;
				$tempa=1;
			}elsif($info[$i] eq "c"){
				$flagb++;
				$tempb=1;
			}elsif($info[$i] eq "d"){
				$flagb--;
				$tempb=1;
			}else{
				print "Error....\n";
			}
		}
		if($flaga>0){
			$tempa = 0;
		}
		if($flagb>0){
			$tempb = 0;
		}
		$DataS{$d} = "$flaga\t$flagb\t$tempa\t$tempb";
	}
	my @OUT;
	my $s=-1;
	my $flag=0;
	foreach my $d (sort{$a<=>$b} keys %DataS){
#               print "Value:$d\t$DataS{$d}\n";
		my @info = split/\t/,$DataS{$d};
		if($out eq "a-b"){
#                       print "a-b $s\t$flag\n";
			if($info[0]>0 && $info[1]==0&& $info[3]==0){
				if($flag==0){
					$flag = 1;
					$s = $d;
				}else{
				}
			}elsif($info[0]>0 && $info[1]==0&& $info[3]==1){
				if ($flag == 1){
					if($s<=$d-1){
						push @OUT,"$s-".($d-1);
					}
				}
				$s = $d+1;
				$flag = 1;

			}elsif($flag==1){
				if($info[1]==0 && $info[3]==0){
					push @OUT,"$s-$d";
				}elsif($s<=$d-1){
					push @OUT,"$s-".($d-1);
				}
				$flag=0;
			}elsif($info[1]==0&&$info[3]==0&&$info[2]>0){
				push @OUT,"$d-$d";
			}
		}elsif($out eq "b-a"){
#                       print "b-a $s\t$flag\n";
			if($info[1]>0 && $info[0]==0&& $info[2]==0){
				if($flag==0){
					$flag = 1;
					$s = $d;
				}else{
				}
			}elsif($info[1]>0 && $info[0]==0&& $info[2]==1){
				if ($flag == 1){
					if($s<=$d-1){
						push @OUT,"$s-".($d-1);
					}
				}
				$s = $d+1;
				$flag = 1;

			}elsif($flag==1){
				if($info[0]==0 && $info[2]==0){
					push @OUT,"$s-$d";
				}elsif($s<=$d-1){
					push @OUT,"$s-".($d-1);
				}
				$flag=0;
			}elsif($info[0]==0&&$info[2]==0&&$info[3]>0){
				push @OUT,"$d-$d";
			}
		}elsif($out eq "aUb"){
#                       print "aUb $s\t$flag\n";
			if($info[0]>0 || $info[1]>0){
				if($flag==0){
					$flag = 1;
					$s = $d;
				}
			}elsif($flag == 1){
				push @OUT,"$s-$d";
				$flag = 0;
			}elsif($flag == 0 &&($info[2]+$info[3])>0){
				push @OUT,"$d-$d";
			}
		}elsif($out eq "aIb"){
#print "aIb $s\t$flag\n";
			if($info[0]>0 && $info[1]>0){
				if($flag==0){
					$flag = 1;
					$s = $d;
				}
			}elsif($flag == 1){
				push @OUT,"$s-$d";
				$flag = 0;
			}elsif($flag == 0 && ($info[0]+$info[2])>0 && ($info[1]+$info[3])>0){
				push @OUT,"$d-$d";
			}

		}else{
			print "Error, Undefined output type $out, should in a-b, b-a, aUb, aIb\n";
		}


	}
	my @OUTorder;
	if(@OUT == 1){
		return \@OUT;
	}else{
		my @info1 = split/-/,$OUT[0];
		my $s = $info1[0];
		for(my$i=1;$i<@OUT;$i++){
			@info1 = split/-/,$OUT[$i-1];
			my @info2 = split/-/,$OUT[$i];
			if($info1[1]+1 == $info2[0]){
				next;
			}else{
				push @OUTorder,"$s-$info1[1]";
				$s = $info2[0];
			}
		}
		@info1 = split/-/,$OUT[-1];
		push @OUTorder,"$s-$info1[1]";
		return \@OUTorder;
	}
}

sub Quantile{
# calculate quantile
# output format: Q0, Q25, Q50, Q75, Q100
	my(@array) = @_;
#	my @array = split/\s+/, $v;
	my @infoS = sort{$a <=> $b} @array;
	my $Q0 = $infoS[0];
	my $Q100 = $infoS[-1];
#	print "$array[0]\t$array[-1]\n";
#	print "$Q0\t$Q100\n";
	my @Q = (25,50,75);
	my %Quantile;
	for(my$i=0;$i<@Q;$i++){
		if ((@infoS % (100/$Q[$i]))!=0) {
			$Quantile{$Q[$i]} = $infoS[int(@infoS/(100/$Q[$i]))];
		}else{
			$Quantile{$Q[$i]} = ($infoS[@infoS/(100/$Q[$i]) - 1] + $infoS[@infoS/(100/$Q[$i])])/2;
		}
	}
	my @QQ = ($Q0, $Quantile{25}, $Quantile{50}, $Quantile{75}, $Q100);
	return \@QQ;
}

sub median {
	my(@array) = @_;
#	my @array = split/\s+/, $m;
	my @infoS = sort{$a <=> $b} @array;
	my $median = 0;
	if ((@infoS % 2)!=0) {
		$median = $infoS[int(@infoS/2)];
	}else{
		$median = ($infoS[@infoS/2 - 1] + $infoS[@infoS/2])/2;
	}
	return $median;
}

sub avg{
	my (@array) = @_;
#	my @info = split/\t/,$str;
	my $sum=0;
	foreach my $i(@array){
		$sum+=$i/@array;
	}
	return $sum;
}

sub sum{
	my (@array) = @_;
	my $sum=0;
	foreach my $i(@array){
		$sum+=$i;
	}
	return $sum;
}

sub min{
	my (@array) = @_;
	my @infoS = sort {$a<=>$b} @array;
	return $infoS[0];
}

sub max{
	my (@array) = @_;
	my @infoS = sort {$a<=>$b} @array;
	return $infoS[-1];
}

sub roundup {
	my $n = shift;
	return(($n == int($n)) ? $n : int($n + 1))
}

1;
