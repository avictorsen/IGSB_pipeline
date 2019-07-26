#!/usr/bin/perl
use strict;
use File::Basename;

if (@ARGV != 1) {
	print "needs genome as argument\n";
	die;
}

###### Get config information
my(%Config);
$Config{'sp'} = @ARGV[0];
$Config{'date'} = `date +%F`;
$Config{'date'} =~ s/\-//g;
chomp $Config{'date'};
$Config{'fastq_suffix'} = 'fastq';
$Config{'bin_dir'} = './idrpipeline/';

my$out_dir="../";
my$bin_dir=$Config{'bin_dir'};
my$prj_dir;
my$fastq_suffix=$Config{'fastq_suffix'};
my$node = 8;
my $IDR_true=0.02;
my $IDR_pr = 0.01;
my$npeaks;
my$genome_table;
my$genome_bedtools;
my$black_list;
my $mito;

if ($Config{'sp'} eq "dm3") {
	$genome_table = "genome_table.Dmel5.41.txt";
	$genome_bedtools = "/dm3/dm3";
        $black_list = "dm3-blacklist.bed";
	$npeaks=30000;
	$prj_dir="/glusterfs/bionimbus/modENCODE_ChIP-seq/";
	$Config{'reference'} = $Config{'bin_dir'}.$genome_bedtools;
	$mito="M";
}elsif ($Config{'sp'} eq "dm6") {
        $genome_table = "genome_table.Dmel6.txt";
        $genome_bedtools = "/dm6/dm6";
        $black_list = "dm6-blacklist.bed";
        $npeaks=30000;
        $prj_dir="/glusterfs/bionimbus/modENCODE_ChIP-seq/";
        $Config{'reference'} = $Config{'bin_dir'}.$genome_bedtools;
	$mito="M";
}elsif ($Config{'sp'} eq "WS220") {
	$genome_table = "genome_table.worm.ws220.txt";
	$genome_bedtools = "/WS220/WS220";
        $black_list = "WS220-blacklist.bed";
        $npeaks=30000;
        $prj_dir="/glusterfs/bionimbus/modENCODE_ChIP-seq/";
        $Config{'reference'} = $Config{'bin_dir'}.$genome_bedtools;
	$mito="chrM";
}elsif ($Config{'sp'} eq "WS245") {
        $genome_table = "genome_table.worm.ws245.txt";
        $genome_bedtools = "/WS245/WS245";
        $black_list = "WS245-blacklist.bed";
        $npeaks=30000;
        $prj_dir="/glusterfs/bionimbus/modENCODE_ChIP-seq/";
        $Config{'reference'} = $Config{'bin_dir'}.$genome_bedtools;
	$mito="CHROMOSOME_M";
}elsif ($Config{'sp'} eq "hg19") {
        $genome_table = "genome_table.human.hg19.txt";
        $genome_bedtools = "/hg19/hg19";
	$black_list = "hg19-blacklist.bed";
        $npeaks=150000;
        $prj_dir="/glusterfs/bionimbus/Human_ENCODE/";
        $Config{'reference'} = $Config{'bin_dir'}.$genome_bedtools;
	$mito="chrM";
}else{
	print "ERROR! Please specifiy genome for this run: hg19; dm3; dm6; WS220; WS245\n";
	die;
}

# paths for tools on Bionimbus cloud
my$bwapath="/home/ubuntu/TOOLS/bwa-0.7.8/";
my$samtoolspath="/home/ubuntu/TOOLS/samtools-0.1.19/";
my$bedtoolspath="/home/ubuntu/TOOLS/bedtools2/bin/";
#my$macspath="/usr/local/bin/";

print "\n#### configure parameters
date\t$Config{'date'}
bin_dir\t$bin_dir
out_dir\t$out_dir
prj_dir\t$prj_dir
node\t$node
IDR_true_threshold\t$IDR_true
IDR_pr_threshold\t$IDR_pr
num_of_start_peaks\t$npeaks

\n";

###############################################################################
open(IN1, "../table.txt") or die $!;
open(OUT1, ">../download.sh") or die $!;
open(OUT2, ">../run_script.sh") or die $!;
open(OUT3, ">../report.sh") or die $!;

###### Get dataset information from data production table
print "#### factors\n";
my$l=0;
my(%Dataset);
while (<IN1>) {
	chomp;
	my @info = split/\t/;
	#print @info;
	next if ($info[7] eq "");
	$l++;
	
	#trim leading and trailing spaces
	my $str;
	foreach (@info) {
		$str= $_;
		$str =~ s/^\s+//;
		$str =~ s/\s+$//;
		$_ = $str;
	}
# Strain_Tissue_STAGE_Source_AbID_Replicate
	my $factor = "$info[1]";
	if ($info[2] ne "")
	{
	 	$factor = "$factor\_$info[2]";
	}
	if ($info[3] ne "")
	{
		$factor = "$factor\_$info[3]";
	}
	my $source = $info[4];
	my $abid;
	if ($source eq "Antibody" || $source eq "IP") {
		$source = "IP";
		$abid = $info[5];
	}
	if ($source eq "Input" || $source eq "INPUT" || $source eq "DNA") {
		$source = "Input";
		$abid = "NA";
	}
	my $replicate;
	if ($info[6] =~ /[Rr]ep/)
	{
		$replicate = $info[6];
	}
	else
	{
		$replicate = "Rep".$info[6];
	} 
	my $name = "$factor\_$source\_$abid\_$replicate";
	
	#either copy raw file from bionimbus or download from url
	
	my $fn = basename($info[7]);
	#print $fn;
	print OUT1 "echo \"--------------------------------------------------------------------------------\"\n";
	if ($info[7] =~ /http/)
	{
		print OUT1 "wget $info[7] .";
	}
	else	
	{	print OUT1 "rsync -c --progress avictorsen\@192\.170\.228\.37\:$prj_dir/$fn ./";
	}
	print OUT1 qq{
for file in ./$fn; do
        format=`file -bi \$file`;
	case \$format in
		\'application\/gzip\'\*\)
		        check=`zcat \$file | awk \'NR \% 4 \=\= 0\' | head -3000 | egrep -c -m 500 \[L-Za-j\]`;
		        if [ \"\$check\" -gt 100 ];then
                		echo "converting gzipped file to Phred33";
		                zcat \$file | perl ./idrpipeline/Phred64_to_Phred33.pl >> $name.$Config{'fastq_suffix'};
				if [ \$? -gt 0 ];then echo "file failed to convert!"; exit; fi;
			else
                                echo "gzipped file is already Phred33";
				zcat \$file >> $name.$Config{'fastq_suffix'};
		        fi;;
		\*\)
                        check=`cat \$file | awk \'NR \% 4 \=\= 0\' | head -3000 | egrep -c -m 500 \[L-Za-j\]`;
		        if [ \"\$check\" -gt 100 ];then
                		echo "converting file to Phred33";
			        cat \$file | perl ./idrpipeline/Phred64_to_Phred33.pl >> $name.$Config{'fastq_suffix'};
                                if [ \$? -gt 0 ];then echo "file failed to convert!"; exit; fi;
			else
                                echo "file is already Phred33";
				cat \$file >> $name.$Config{'fastq_suffix'};
			fi;;
	esac;
	if [ -e $name.$Config{'fastq_suffix'}.gz ];then
		zcat $name.$Config{'fastq_suffix'}.gz >> $name.$Config{'fastq_suffix'};
		if [ \$? -gt 0 ];then exit; else rm $name.$Config{'fastq_suffix'}.gz; fi;
		pigz $name.$Config{'fastq_suffix'} || { exit ;};
	else
		pigz $name.$Config{'fastq_suffix'} || { exit ;};
	fi;
	rm \$file;
done\n\n};
	if (!defined $Dataset{$factor}{$source}{$replicate}) {
	#	print OUT1 "\nmv $fn $name.$Config{'fastq_suffix'} || { exit ;}\n\n";
		$Dataset{$factor}{$source}{$replicate}=$name;
	;
	}else{
	#	print OUT1 "\ncat $fn >> $Dataset{$factor}{$source}{$replicate}.$Config{'fastq_suffix'}\n";
	#	print OUT1 "rm -f $fn\n\n";
		print "!!!!!!!!!@@@@@@@@ WARNING: $factor\_$source\_$abid\_Rep$replicate has duplicated files, will be merged for processing\n";
	}
	print "$factor\t$source\t$abid\t$replicate\t$fn\n";
}
print OUT1 "mail -s 'download finished on '\$HOSTNAME -aFROM:\\<avictorsen\@IGSB.opensciencedatacloud.org\\> avictorsen\@gmail.com <<< 'this message has no content.'";
close IN1;

my$dn=keys %Dataset;
print "\n\n#### Summary: 
$dn Datasets were imported!
$l fastq files are being processed!
\n";

###### Organize data sets, do alignment
my(%Peakset);
print OUT2 "############ Organize data sets, and run alignment\n";
#print OUT2 "ln -fs $Config{'bin_dir'}/genome_tables/$genome_table genome_table.txt\n";
#print OUT2 "ln -fs $Config{'bin_dir'}/*.r ./\n";
print OUT2 "cp $Config{'bin_dir'}/genome_tables/$genome_table genome_table.txt\n";
print OUT2 "cp $Config{'bin_dir'}/*.r ./\n";
foreach my $d (sort {$a cmp $b} keys %Dataset) {
	foreach my $r (sort {$a <=> $b} keys %{$Dataset{$d}{'IP'}}) {
		my$ip = $Dataset{$d}{'IP'}{$r};
		$Peakset{$d}.="$ip\t";
	}
	print OUT2 "\n###### Dataset $d
	echo \"MSG: running alignment for $d\"
	#gzip -d $d\*
	## aligning $d and write sam file
	#for file in $d\*\.$fastq_suffix; do
	#	$bwapath/bwa aln -B 0 -t $node $Config{'reference'} \$file | $bwapath/bwa samse $Config{'reference'} - \$file > \$\{file\%fastq}sam\;
        #        if [ \$? -gt 0 ];then exit; else rm \$file; fi;
	#done
	## create bamfiles
		find \. \-name \"$d\*\.bam\" | xargs \-n 1 \-P $node \-iFILES bash -c \'
			out=\`echo FILES | sed -e \"s\/bam\$\/sam\/g\"\`
			$samtoolspath/samtools view -h FILES \> \$out
		\'\;
";
	## if genome is WS245 convert sam file chromosomes and remake BAM files
	if ($Config{'sp'} eq 'WS245'){
		print OUT2 "\n\t## rewrite chromosomes and remake bam files
		for file in $d\*\.sam; do
			awk \-F\"\\t\" \'\{if \(\$0 \~ \/\^\@\/\)\{sub\(\/SN\\\:chr\/,\"SN\:CHROMOSOME\_\",\$0\);print \$0;\} else\{if \(\$3 \~ \/chr\/\)\{sub\(\/chr\/,\"CHROMOSOME_\",\$3\)\} else\{\;\}\;print \$0\;\}\}\' OFS\=\"\\t\" \$file \> \$\{file\%sam\}temp.sam
			mv \$\{file\%sam\}temp.sam \$file
		done
		find \. \-name \"$d\*\.sam\" | xargs \-n 1 \-P $node \-iFILES bash -c \'
			out=\`echo FILES | sed -e \"s\/sam\$\/bam\/g\"\`
			$samtoolspath/samtools view -Sbh FILES \> \$out
		\'\;
";
	}
	print OUT2 "\r
        ## create q30 file
		find \. \-name \"$d\*\.sam\" | xargs \-n 1 \-P $node \-iFILES bash -c \'
			out=\`echo FILES | sed -e \"s\/\.sam\$\/\_q30\.sam\/g\"\`
			$samtoolspath/samtools view -F 1548 -h -S -q 30 -o \$out FILES;
		\'\;
	## create tagAlign files
		find \. \-name \"$d\*q30\.sam\" | xargs \-n 1 \-P $node \-iFILES bash -c \'
			out=\`echo FILES | sed -e \"s\/\_q30\.sam\$\/\.tagAlign\.gz\/g\"\`
			$samtoolspath/samtools view -Sb FILES | $bedtoolspath/bamToBed -i stdin | awk \'\\\'\'BEGIN\{FS\=\"\\t\"\;OFS=\"\\t\"\}\{\$4\=\"N\"\;print\$0\}\'\\\'\' | gzip -c \> \$out
		\'\;
	## create pooled IPs and INPUTs
		for file in ./$d\_IP\_\*\.tagAlign\.gz; do
			zcat \$file >> $d\_IP_Rep0.tagAlign;
		done
		pigz $d\_IP_Rep0.tagAlign

		for file in ./$d\_Input\_\*\.tagAlign\.gz; do
			zcat \$file >> $d\_Input_Rep0.tagAlign;
		done
		pigz $d\_Input_Rep0.tagAlign
	## create pseudoreps for IPs
		find \. \-name \"$d\*\_IP\_\*tagAlign\.gz\" | xargs \-n 1 \-P $node \-iFILES bash -c \'
			out=`echo FILES | sed -e \"s\/\.tagAlign\.gz\/\/g\"`
			./idrpipeline/shuf.sh FILES ./ \$out
		\'\;
	## generate pbc quality control files
		find \. \-name \"$d\*\.bam\" | xargs \-n 1 \-P $node \-iFILES bash -c \'
			$bedtoolspath/bedtools bamtobed -i FILES | awk \'\\\'\'BEGIN{OFS=\"\\t\"}{print \$1,\$2,\$3,\$6}\'\\\'\' | grep -v $mito | sort | uniq -c | awk \'\\\'\'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1\+1} (\$1==2){m2=m2\+1} {m0=m0\+1} {mt=mt\+\$1} END{printf \"\%d\\t\%d\\t\%d\\t\%d\\t\%f\\t\%f\\t\%f\\n\",mt,m0,m1,m2,m0\/mt,m1\/m0,m1\/m2}\'\\\'\' > FILES.pbc.qc
		\'\;
	## generate Rep0.bam
		$bedtoolspath/bedToBam -i $d\_Input_Rep0.tagAlign.gz -g $Config{'bin_dir'}/$genome_bedtools.genome > $d\_Input_Rep0.bam
	## generate SAM stats
		perl $Config{'bin_dir'}/SAM.statistics.bwa.pl ./ $Config{'date'}.SAM.stats
";

}
print OUT2 "\n########### remove sam files
for file in \.\/\*\.sam;do
	if [ -e \$\{file\%sam\}bam ] ; then
		rm -f \$file;
	fi;
done";

###### run SPP and IDR
print OUT2 "\n############ SPP and IDR\n";
foreach my $d (sort {$a cmp $b} keys %Peakset) {
	print OUT2 "###### Dataset\t$d
	echo \"MSG: running SPP and IDR for $d\"
	echo \"MSG: IPs: $Peakset{$d}\"
	echo \"MSG: Inputs: $d\_Input_Rep0.tagAlign.gz\"
	";
	print OUT3 "###### Dataset\t$d
	echo \"IDR-filtered SPP/MACS2 $d\"
	echo \"IPs:\t$Peakset{$d}\"
	echo \"Inputs:\t$d\_Input_Rep0.tagAlign.gz\"
	max_numPeaks_Rep=0
	#max_numPeaks_Rep_macs2=0
";
	my@temp=split/\t/,$Peakset{$d};
	my@IPs=sort { $a cmp $b } @temp;
	my(@IDR_pairs,@IDR_pairs_macs);
	for(my$i=0;$i<@IPs;$i++){
		print OUT2 "# peak calls: original replicates and self-pesduoreplicates $IPs[$i]
		Rscript $Config{'bin_dir'}/run_spp.wig.R -c=$IPs[$i].tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc
		pigz -cdf $IPs[$i].tagAlign.gz > $IPs[$i].tagAlign
		Rscript $bin_dir/spp_wig_tagAlign.R  --args  --baseDir=./ --outDir=./ --factor=\"$IPs[$i].tagAlign\" --ipFile=$IPs[$i].tagAlign
		Rscript $bin_dir/run_spp.R -c=$IPs[$i].pr1.tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc
		Rscript $Config{'bin_dir'}/run_spp.R -c=$IPs[$i].pr2.tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc
	";
	}
	print OUT2 "# peak calls on pooled rep
		Rscript $Config{'bin_dir'}/run_spp.wig.R -c=$d\_IP_Rep0.tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc
		pigz -cdf $d\_IP_Rep0.tagAlign.gz > $d\_IP_Rep0.tagAlign
		Rscript $bin_dir/spp_wig_tagAlign.R  --args  --baseDir=./ --outDir=./ --factor=\"$d\_IP_Rep0.tagAlign\" --ipFile=$d\_IP_Rep0.tagAlign
	# peak calls on pooled pseudo replicates
		Rscript $Config{'bin_dir'}/run_spp.R -c=$d\_IP_Rep0.pr1.tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc
		Rscript $Config{'bin_dir'}/run_spp.R -c=$d\_IP_Rep0.pr2.tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc
	# peak calls on input
		Rscript $bin_dir/spp_wig_BAM.R  --args  --baseDir=./ --outDir=./ --factor=$d\_Input_Rep0 --ipFile=$d\_Input_Rep0.bam
	#####################
		perl $Config{'bin_dir'}/remove_1strow_encodepeak.pl ./ ./\n";

	## consistancy
	if (scalar @IPs <= $node){
		my$combined=substr($IPs[0],0,-1);
		print OUT2 "\t# consistency: self and pesduoreplicates
		IPs=\"$d\_IP_Rep0.pr1.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz;";
		for (my$i=0;$i<@IPs;$i++){
			print OUT2 "$IPs[$i].pr1.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz;";
		}
		print OUT2 "\"
		echo -n \$IPs | xargs -d\"\;\" -P 8 -iFILES bash -c '
			FILE1=`echo FILES`
			FILE2=`echo FILES | sed \'s\/\.pr1\.tagAlign\/\.pr2\.tagAlign\/\'`
			OUTPUT=\${FILE1%.pr1.tagAlign_VS_*}
			Rscript batch-consistency-analysis.r \$FILE1 \$FILE2 -1 \$OUTPUT.pr1_VS_\$OUTPUT.pr2 0 F signal.value
			Rscript batch-consistency-plot.r 1 \$OUTPUT.pr1_VS_\$OUTPUT.pr2 \$OUTPUT.pr1_VS_\$OUTPUT.pr2
		\'\;";
                for (my$i=0;$i<@IPs;$i++){
		}
		print OUT2 "\n\t# consistency: original replicates
		IPs=\"";
		for(my$i=0;$i<@IPs;$i++){
			for(my$j=$i+1;$j<@IPs;$j++){
				print OUT2 "$IPs[$j].tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz $IPs[$i].tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz $IPs[$j]_VS_$IPs[$i];";
			}
		}
		print OUT2 "\"
		echo -n \$IPs | xargs -d\"\;\" -P8 -iFILES bash -c '
			FILE1=`echo FILES | cut -d\" \" -f1`
			FILE2=`echo FILES | cut -d\" \" -f2`
			OUTPUT=`echo FILES | cut -d\" \" -f3`
			Rscript batch-consistency-analysis.r \$FILE1 \$FILE2 -1 \$OUTPUT 0 F signal.value
		\'\;
";
		print OUT2 "\t# Generate bed files
		max_numPeaks_Rep=0
		#max_numPeaks_Rep_macs2=0
";
                for (my$i=0;$i<@IPs;$i++){
                        my $R1=substr($IPs[$i],-1);
                        print OUT2 "\t\t## IDR with Rep$R1
";
                        print OUT2 "\t\t\tnumPeaks_Rep$R1\_pr=\$(awk '\$11 <= $IDR_pr {print \$0}' $IPs[$i].pr1_VS_$IPs[$i].pr2-overlapped-peaks.txt | wc -l)\n";
                        for (my$j=$i+1;$j<@IPs;$j++){
                                my$R2=substr($IPs[$j],-1);
                                print OUT2 "\n\t\t\tnumPeaks_Rep$R2\_Rep$R1=\$(awk '\$11 <= $IDR_true {print \$0}' $IPs[$j]_VS_$IPs[$i]-overlapped-peaks.txt | wc -l)\n";
                                print OUT2 "\t\t\tif [ \$numPeaks_Rep$R2\_Rep$R1 -gt \$max_numPeaks_Rep ]; then\n\t\t\t\tmax_numPeaks_Rep=\$numPeaks_Rep$R2\_Rep$R1\n\t\t\tfi\n";
                        }
                }
                print OUT2 "\t\t## IDR for combined Rep
";
                print OUT2 "\t\t\tnumPeaks_Rep0_pr=\$(awk '\$11 <= $IDR_pr {print \$0}' $d\_IP_Rep0.pr1_VS_$d\_IP_Rep0.pr2-overlapped-peaks.txt | wc -l)\n";
                print OUT2 "\t\t## print optimal and conservitive IDR values
";
                print OUT2 "\t\t\tif [ \$numPeaks_Rep0_pr -gt \$max_numPeaks_Rep ]; then\n\t\t\t\toptThresh=\$numPeaks_Rep0_pr\n\t\t\telse\n\t\t\t\toptThresh=\$max_numPeaks_Rep\n\t\t\tfi\n";
		print OUT2 "\t\tzcat $d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz | sort -k7nr,7nr | head -n \$optThresh | pigz -c > spp.optimal.$d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz\n";
		print OUT2 "\t\tzcat $d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz | sort -k7nr,7nr | head -n \$max_numPeaks_Rep | pigz -c >  spp.conservative.$d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz\n";

		# print lines for report.sh
		for (my$i=0;$i<@IPs;$i++){
			my $R1=substr($IPs[$i],-1);
			print OUT3 "\t## IDR with Rep$R1
";
			print OUT3 "\t\tnumPeaks_Rep$R1\_pr=\$(awk '\$11 <= $IDR_pr {print \$0}' $IPs[$i].pr1_VS_$IPs[$i].pr2-overlapped-peaks.txt | wc -l)\n\t\techo \"numPeaks_Rep$R1\_pr\t\$numPeaks_Rep$R1\_pr\"\n";
			for (my$j=$i+1;$j<@IPs;$j++){
				my$R2=substr($IPs[$j],-1);
				print OUT3 "\n\t\tnumPeaks_Rep$R2\_Rep$R1=\$(awk '\$11 <= $IDR_true {print \$0}' $IPs[$j]_VS_$IPs[$i]-overlapped-peaks.txt | wc -l)\n\t\techo \"numPeaks_Rep$R2\_Rep$R1\t\$numPeaks_Rep$R2\_Rep$R1\"\n";
				print OUT3 "\t\tif [ \$numPeaks_Rep$R2\_Rep$R1 -gt \$max_numPeaks_Rep ]; then\n\t\t\tmax_numPeaks_Rep=\$numPeaks_Rep$R2\_Rep$R1\n\t\tfi\n";
			}
		}
		print OUT3 "\t## IDR for combined Rep
";
		print OUT3 "\t\tnumPeaks_Rep0_pr=\$(awk '\$11 <= $IDR_pr {print \$0}' $d\_IP_Rep0.pr1_VS_$d\_IP_Rep0.pr2-overlapped-peaks.txt | wc -l)\n\t\techo \"numPeaks_Rep0_pr\t\$numPeaks_Rep0_pr\"\n";
		print OUT3 "\t## print optimal and conservitive IDR values
";
		print OUT3 "\t\tif [ \$numPeaks_Rep0_pr -gt \$max_numPeaks_Rep ]; then\n\t\t\toptThresh=\$numPeaks_Rep0_pr\n\t\telse\n\t\t\toptThresh=\$max_numPeaks_Rep\n\t\tfi\n\t\techo \"optThresh\t\$optThresh\"\n\t\techo \"conThresh\t\$max_numPeaks_Rep\"\n";
	}
	else{die "too many IP replicates";}
}

print OUT2 "
########## Organizing files
### move wiggle files to subdirectory
mkdir wig/
mv *.wig wig/

### generate IDR report
sh report.sh > $Config{'sp'}_idr.log

### generate bed files for SPP peak calls
mkdir bed/
ls *regionPeak.gz > regionPeak.list
while read -r line; do
	zcat \$line | cut -f1-3 > \$line.bed33
done < regionPeak.list
mv *regionPeak.gz.bed33 bed/
rename s/regionPeak.gz.bed33/bed/ bed/*

### removed unzipped files
rm -f *tagAlign

### remove 50% overlapping peaks from blacklist
ls spp.*.regionPeak.gz  > spp.list
while read -r line; do
	bedtools intersect -a \$line -b $bin_dir/$black_list  -f 0.5 -v > \$line.rmblacklist
done < spp.list
rename s/regionPeak.gz.rmblacklist/regionPeak.rmblacklist/ *.regionPeak.gz.rmblacklist

### write final files to STDOUT to capture in nohup.txt
cat $Config{'date'}.SAM.stats
cat $Config{'date'}.cc
cat $Config{'sp'}_idr.log

### report as finished via email
mail -s 'Run on '\$HOSTNAME' has finished!' -aFROM:\\<avictorsen\@IGSB.opensciencedatacloud.org\\> avictorsen\@gmail.com < $Config{'sp'}_idr.log
";

close OUT1;
close OUT2;
close OUT3;
