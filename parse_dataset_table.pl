#!/usr/bin/perl
use strict;
use File::Basename;

if (@ARGV != 3) {
	print "This program takes arguments:\n\trun bwa: Y/N\n\tgenome: dm3, dm6, hg19, WS220, WS245\n\tpeak caller: macs2, spp, or both\n";
	exit;
}

###### Get config information
my(%Config);
$Config{'aln'} = @ARGV[0];
$Config{'sp'} = @ARGV[1];
$Config{'peak'} = @ARGV[2];
$Config{'date'} = `date +%F`;
$Config{'date'} =~ s/\-//g;
chomp $Config{'date'};
$Config{'fastq_suffix'} = 'fastq';
$Config{'bin_dir'} = '../sf_idrpipeline/';

my$out_dir="./";
my$bin_dir=$Config{'bin_dir'};
my$prj_dir;
my$fastq_suffix=$Config{'fastq_suffix'};
my$node = 8;
my$IDR_true=0.02;
my$IDR_pr = 0.01;
my$npeaks;
my$genome_table;
my$genome_bedtools;
my$black_list;
my$mito;
my$initials;

if ($Config{'peak'} ne "macs2" && $Config{'peak'} ne "spp" && $Config{'peak'} ne "both"){
	print "3nd argument needs to be:\n\tmacs2: to use macs2 peak caller\n\tspp: to use spp peak caller\n\tboth: to use both macs2 and spp peak callers\n";
	exit;
}

if ($Config{'sp'} eq "dm3") {
	$genome_table = "genome_table.Dmel5.41.txt";
	$genome_bedtools = "/dm3/dm3";
        $black_list = "dm3-blacklist.bed";
	$npeaks=30000;
	#$prj_dir="/glusterfs/bionimbus/modENCODE_ChIP-seq/";
	$Config{'reference'} = $Config{'bin_dir'}.$genome_bedtools;
	$mito="M";
	$initials="dm";
}elsif ($Config{'sp'} eq "dm6" || $Config{'sp'} eq "dm6.masked") {
        $genome_table = "genome_table.Dmel6.txt";
	if ($Config{'sp'} eq "dm6.masked"){
	        $genome_bedtools = "/dm6.masked/dm6.masked";
	}else{
	        $genome_bedtools = "/dm6/dm6";
	}
        $black_list = "dm6-blacklist.bed";
        $npeaks=30000;
        #$prj_dir="/glusterfs/bionimbus/modENCODE_ChIP-seq/";
        $Config{'reference'} = $Config{'bin_dir'}.$genome_bedtools;
	$mito="M";
	$initials="dm";
}elsif ($Config{'sp'} eq "WS220") {
	$genome_table = "genome_table.worm.ws220.txt";
	$genome_bedtools = "/WS220/WS220";
        $black_list = "WS220-blacklist.bed";
        $npeaks=30000;
        #$prj_dir="/glusterfs/bionimbus/modENCODE_ChIP-seq/";
        $Config{'reference'} = $Config{'bin_dir'}.$genome_bedtools;
	$mito="chrM";
	$initials="ce";
}elsif ($Config{'sp'} eq "WS245") {
        $genome_table = "genome_table.worm.ws245.txt";
        $genome_bedtools = "/WS245/WS245";
        $black_list = "WS245-blacklist.bed";
        $npeaks=30000;
        #$prj_dir="/glusterfs/bionimbus/modENCODE_ChIP-seq/";
        $Config{'reference'} = $Config{'bin_dir'}.$genome_bedtools;
	$mito="CHROMOSOME_M";
	$initials="ce";
}elsif ($Config{'sp'} eq "hg19") {
        $genome_table = "genome_table.human.hg19.txt";
        $genome_bedtools = "/hg19/hg19";
	$black_list = "hg19-blacklist.bed";
        $npeaks=150000;
        #$prj_dir="/glusterfs/bionimbus/Human_ENCODE/";
        $Config{'reference'} = $Config{'bin_dir'}.$genome_bedtools;
	$mito="chrM";
	$initials="ho";
}else{
	print "ERROR! Please specifiy genome as second argument for this run: hg19; dm3; dm6; WS220; WS245\n";
	exit;
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
node\t$node
IDR_true_threshold\t$IDR_true
IDR_pr_threshold\t$IDR_pr
num_of_start_peaks\t$npeaks

\n";

###############################################################################
open(IN1, "./table.txt") or die $!;
open(OUT1, ">./download.sh") or die $!;
open(OUT2, ">./run_script.sh") or die $!;
open(OUT3, ">./report.sh") or die $!;

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
#	if ($info[7] =~ /http/)
#	{
#		print OUT1 "wget $info[7] .";
#	}
#	else	
#	{	print OUT1 "rsync -c --progress avictorsen\@192\.170\.228\.37\:$prj_dir/$fn ./";
#	}
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
print OUT2 "cp $Config{'bin_dir'}/genome_tables/$genome_table genome_table.txt\n";
print OUT2 "cp $Config{'bin_dir'}/*.r ./\n";
foreach my $d (sort {$a cmp $b} keys %Dataset) {
	foreach my $r (sort {$a <=> $b} keys %{$Dataset{$d}{'IP'}}) {
		my$ip = $Dataset{$d}{'IP'}{$r};
		$Peakset{$d}.="$ip\t";
	}
	print OUT2 "\n###### Dataset $d\n";
	if ($Config{'aln'} eq "Y"){
		printf OUT2 ("\t%s\n\t%s\n\t%s\n\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t%s\n\t%s\n\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t%s\n\t%s\n\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t%s\n",
			"echo \"MSG: running alignment for $d\"",
			"gzip -d $d\*",
			"## aligning $d and write sam file",
			"for file in $d\*\.$fastq_suffix; do",
				"$bwapath/bwa aln -B 0 -t $node $Config{'reference'} \$file | $bwapath/bwa samse $Config{'reference'} - \$file > \$\{file\%fastq}sam\;",
                		"if [ \$? -gt 0 ];then exit; else rm \$file; fi;",
			"done",
			"## create q30 file",
			"find \. \-name \"$d\*\.sam\" | xargs \-n 1 \-P $node \-iFILES bash -c \'",
				"out=\`echo FILES | sed -e \"s\/\.sam\$\/\_q30\.sam\/g\"\`",
				"$samtoolspath/samtools view -F 1548 -h -S -q 30 -o \$out FILES;",
			"\'\;",
			"## create bamfiles",
			"find \. \-name \"$d\*\.sam\" | xargs \-n 1 \-P $node \-iFILES bash -c \'",
				"out=\`echo FILES | sed -e \"s\/sam\$\/bam\/g\"\`",
				"$samtoolspath/samtools view -Sbh FILES \> \$out",
			"\'\;"
		);
	}
	elsif ($Config{'aln'} eq "N"){
		printf OUT2 ("\t%s\n\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t%s\n\t%s\n\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t%s\n",
			"## recreate q30.bam files",
			"find \. \-name \"$d\*\.bam\" | xargs \-n 1 \-P $node \-iFILES bash -c \'",
				"out=\`echo FILES | sed -e \"s\/\.bam\$\/\_q30\.bam\/g\"\`",
				"$samtoolspath/samtools view -F 1548 -h -q 30 -b -o \$out FILES",
			"\'\;",
			"## recreate sam files",
			"find \. \-name \"$d\*\.bam\" | xargs \-n 1 \-P $node \-iFILES bash -c \'",
				"out=\`echo FILES | sed -e \"s\/bam\$\/sam\/g\"\`",
				"$samtoolspath/samtools view -h FILES \> \$out",
			"\'\;"
		);
		## if genome is WS245 convert sam file chromosomes and remake BAM files",
		if ($Config{'sp'} eq 'WS245' | $Config{'sp'} eq 'dm3'){
			printf OUT2 ("\t%s\n\t\t%s\n",
				"## rewrite chromosomes and remake bam files",
				"for file in $d\*\.sam; do",
			);
			if ($Config{'sp'} eq 'WS245'){
				printf OUT2 ("\t\t\t%s\n","awk \-F\"\\t\" \'\{if \(\$0 \~ \/\^\@\/\)\{sub\(\/SN\\\:chr\/,\"SN\:CHROMOSOME\_\",\$0\);print \$0;\} else\{if \(\$3 \~ \/chr\/\)\{sub\(\/chr\/,\"CHROMOSOME_\",\$3\)\} else\{\;\}\;print \$0\;\}\}\' OFS\=\"\\t\" \$file \> \$\{file\%sam\}temp.sam");
			}
			if ($Config{'sp'} eq 'dm3'){
				printf OUT2 ("\t\t\t%s\n","awk \-F\"\\t\" \'\{if \(\$0 \~ \/\^\@\/\)\{sub\(\/SN\\\:chr\/,\"SN\:\",\$0\);print \$0;\} else\{if \(\$3 \~ \/chr\/\)\{sub\(\/chr\/,\"\",\$3\)\} else\{\;\}\;print \$0\;\}\}\' OFS\=\"\\t\" \$file \> \$\{file\%sam\}temp.sam");
			}
			printf OUT2 ("\t\t\t%s\n\t\t%s\n\t%s\n\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t%s\n",
					"mv \$\{file\%sam\}temp.sam \$file",
				"done",
				"## rewrite corrected bam file",
				"find \. \-name \"$d\*\.sam\" | xargs \-n 1 \-P $node \-iFILES bash -c \'",
				"out=\`echo FILES | sed -e \"s\/sam\$\/bam\/g\"\`",
				"$samtoolspath/samtools view -Sbh FILES \> \$out",
			"\'\;",
			);
		}
	}
	print OUT2 "\t## generate pbc quality control files and create tagAlign files
		find \. \-name \"$d\*\.bam\" | xargs \-n 1 \-P $node \-iFILES bash -c \'
			OUT=`echo FILES | sed -e \"s\/\_q30\.bam\$\/\.tagAlign\.gz\/g\"`
			$bedtoolspath/bedtools bamtobed -i FILES > FILES.bed
			awk \'\\\'\'BEGIN{OFS=\"\\t\"}{print \$1,\$2,\$3,\$6}\'\\\'\' FILES.bed | grep -v $mito | sort | uniq -c | awk \'\\\'\'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1\+1} (\$1==2){m2=m2\+1} {m0=m0\+1} {mt=mt\+\$1} END{printf \"\%d\\t\%d\\t\%d\\t\%d\\t\%f\\t\%f\\t\%f\\n\",mt,m0,m1,m2,m0\/mt,m1\/m0,m1\/m2}\'\\\'\' > FILES.pbc.qc
			if [[ FILES = *_q30.bam ]]; then
				awk \'\\\'\'BEGIN\{FS\=\"\\t\"\;OFS=\"\\t\"\}\{\$4\=\"N\"\;print\$0\}\'\\\'\' FILES.bed | gzip -c \> \$OUT
			fi
			rm -f FILES.bed
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
			$Config{'bin_dir'}/shuf.sh FILES ./ \$out
		\'\;
	## generate Rep0.bam
		$bedtoolspath/bedToBam -i $d\_Input_Rep0.tagAlign.gz -g $Config{'bin_dir'}/$genome_bedtools.genome > $d\_Input_Rep0.bam
	## generate SAM stats
		perl $Config{'bin_dir'}/SAM.statistics.bwa.pl ./ $Config{'date'}.SAM.stats
";

}

#already removed
#print OUT2 "\n########### remove sam files
#for file in \.\/\*\.sam;do
#	if [ -e \$\{file\%sam\}bam ] ; then
#		rm -f \$file;
#	fi;
#done";

###### run SPP and IDR
print OUT2 "\n############ Call Peaks and IDR\n";
foreach my $d (sort {$a cmp $b} keys %Peakset) {
	printf OUT2 ("%s\n%s\n%s\n\%s\n",
		"###### Dataset\t$d",
		"\techo \"MSG: calling peaks and IDR for $d\"",
		"\techo \"MSG: IPs: $Peakset{$d}\"",
		"\techo \"MSG: Inputs: $d\_Input_Rep0.tagAlign.gz\""
	);
	printf OUT3 ("%s\n%s\n%s\n\%s\n",
		"###### Dataset\t$d",
		"\techo \"IDR-filtered SPP/MACS2 $d\"",
		"\techo \"IPs:\t$Peakset{$d}\"",
		"\techo \"Inputs:\t$d\_Input_Rep0.tagAlign.gz\""
	);
	my@temp=split/\t/,$Peakset{$d};
	my@IPs=sort { $a cmp $b } @temp;
	my(@IDR_pairs,@IDR_pairs_macs);
	printf OUT2 ("\t%s\n\t\t%s\n",
		"# generate input wig",
		"Rscript $bin_dir/spp_wig_BAM.R  --args  --baseDir=./ --outDir=./ --factor=$d\_Input_Rep0 --ipFile=$d\_Input_Rep0.bam"
	);

	if ($Config{'peak'} eq "spp" || $Config{'peak'} eq "both"){
		for(my$i=0;$i<@IPs;$i++){
			printf OUT2 ("\t%s\n\t\t%s\n\t\t%s\n\t\t%s\n\t\t%s\n\t\t%s\n",
				"# spp peak calls: original replicates and self-pesduoreplicates $IPs[$i]",
				"Rscript $Config{'bin_dir'}/run_spp.wig.R -c=$IPs[$i].tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc",
				"pigz -cdf $IPs[$i].tagAlign.gz > $IPs[$i].tagAlign",
				"Rscript $bin_dir/spp_wig_tagAlign.R  --args  --baseDir=./ --outDir=./ --factor=\"$IPs[$i].tagAlign\" --ipFile=$IPs[$i].tagAlign",
				"Rscript $bin_dir/run_spp.R -c=$IPs[$i].pr1.tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc",
				"Rscript $Config{'bin_dir'}/run_spp.R -c=$IPs[$i].pr2.tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc"
			);
		}
			printf OUT2 ("\t%s\n\t\t%s\n\t\t%s\n\t\t%s\n\t%s\n\t\t%s\n\t\t%s\n",
				"# spp peak calls on pooled rep",
				"Rscript $Config{'bin_dir'}/run_spp.wig.R -c=$d\_IP_Rep0.tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc",
				"pigz -cdf $d\_IP_Rep0.tagAlign.gz > $d\_IP_Rep0.tagAlign",
				"Rscript $bin_dir/spp_wig_tagAlign.R  --args  --baseDir=./ --outDir=./ --factor=\"$d\_IP_Rep0.tagAlign\" --ipFile=$d\_IP_Rep0.tagAlign",
				"# spp peak calls on pooled pseudo replicates",
				"Rscript $Config{'bin_dir'}/run_spp.R -c=$d\_IP_Rep0.pr1.tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc",
				"Rscript $Config{'bin_dir'}/run_spp.R -c=$d\_IP_Rep0.pr2.tagAlign.gz -i=$d\_Input_Rep0.tagAlign.gz -npeak=$npeaks -x=-500:85 -odir=./ -p=$node -filtchr='.*_[CD].*' -savr -savp -rf -out=$Config{'date'}.cc"
			);
	}
	if ($Config{'peak'} eq "macs2" || $Config{'peak'} eq "both"){
		printf OUT2 ("\t%s\n\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t%s\n",
			"# macs2 peak calls: original replicates and self-pesduoreplicates $d",
			"ls \-1 $d\_IP\*\.tagAlign\.gz | xargs \-n 1 \-P $node \-iFILES bash -c \'",
				"INPUT=`echo $d\_Input\_Rep0\.tagAlign\.gz`",
				"OUT=`echo FILES | sed -e \"s\/\.tagAlign\.gz\/\.tagAlign\\_VS\\_\$\{INPUT\%\\.gz}\/g\"`",
				"ss=`expr \$(cat $Config{'date'}.cc | grep FILES | cut -f3 | head -n 1 | cut -d \",\" -f1) \/ 2`",
				"macs2 callpeak -t FILES -c \$INPUT -g $initials --to-large --nomodel --extsize \$ss -p 0.1 -n \$OUT",
				"sort -k8nr,8nr \$\{OUT\}\\_peaks\.narrowPeak | head -n $npeaks | gzip -c \> \$\{OUT%\_peaks.narrowPeak}\.narrowPeak\.gz",
				"if [ -e \$\{OUT%\_peaks.narrowPeak}\.narrowPeak\.gz ]; then rm -f \$\{OUT\}\\_peaks\.narrowPeak; fi",
			"\'\;",
		);
	}
	printf OUT2 ("\t%s\n\t\t%s\n\t%s\n",
		"#####################",
		"perl $Config{'bin_dir'}/remove_1strow_encodepeak.pl ./ ./",
		"#####################"
	);
	## consistancy
	if (scalar @IPs <= $node){
		my$combined=substr($IPs[0],0,-1);
		print OUT2 "\t# consistency: self and pesduoreplicates
		IPs=\"$d\_IP_Rep0.pr1.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz;";
		for (my$i=0;$i<@IPs;$i++){
			print OUT2 "$IPs[$i].pr1.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz;";
		}
		printf OUT2 ("\"\n\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n",
			"echo -n \$IPs | xargs -d\"\;\" -P 8 -iFILES bash -c '",
			"FILE1=`echo FILES`",
			"FILE2=`echo FILES | sed \'s\/\.pr1\.tagAlign\/\.pr2\.tagAlign\/\'`",
			"OUTPUT=\${FILE1%.pr1.tagAlign_VS_*}"
		);
		if ($Config{'peak'} eq "spp" || $Config{'peak'} eq "both"){
			printf OUT2 ("\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n",
				"#spp",
				"Rscript batch-consistency-analysis.r \$FILE1 \$FILE2 -1 \$OUTPUT.pr1_VS_\$OUTPUT.pr2 0 F signal.value",
				"Rscript batch-consistency-plot.r 1 \$OUTPUT.pr1_VS_\$OUTPUT.pr2 \$OUTPUT.pr1_VS_\$OUTPUT.pr2"
			);
		}
		if ($Config{'peak'} eq "macs2" || $Config{'peak'} eq "both"){
			printf OUT2 ("\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n",
				"#macs2",
				"FILE1macs=`echo FILES | sed \'s\/\.regionPeak\.gz\/\.narrowPeak\.gz\/\'`",
				"FILE2macs=`echo FILES | sed \'s\/\.pr1\.tagAlign\/\.pr2\.tagAlign\/\' | sed \'s\/\.regionPeak\.gz\/\.narrowPeak\.gz\/\'`",
				"Rscript batch-consistency-analysis.r \$FILE1macs \$FILE2macs -1 \$OUTPUT.pr1_VS_\$OUTPUT.pr2.macs 0 F p.value",
				"Rscript batch-consistency-plot.r 1 \$OUTPUT.pr1_VS_\$OUTPUT.pr2.macs \$OUTPUT.pr1_VS_\$OUTPUT.pr2.macs"
			);
		}
		print OUT2 "\t\t\'\;";
                for (my$i=0;$i<@IPs;$i++){
		}
		print OUT2 "\n\t# consistency: original replicates
		IPs=\"";
		for(my$i=0;$i<@IPs;$i++){
			for(my$j=$i+1;$j<@IPs;$j++){
				print OUT2 "$IPs[$j].tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz $IPs[$i].tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz $IPs[$j]_VS_$IPs[$i];";
			}
		}
		print OUT2 "\"";

		printf OUT2 ("\n\t\t%s\n\t\t\t%s\n",
			"echo -n \$IPs | xargs -d\"\;\" -P8 -iFILES bash -c '",
				"OUTPUT=`echo FILES | cut -d\" \" -f3`"
		);
		if ($Config{'peak'} eq "spp" || $Config{'peak'} eq "both"){
			printf OUT2 ("\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n",
				"#spp",
				"FILE1=`echo FILES | cut -d\" \" -f1`",
				"FILE2=`echo FILES | cut -d\" \" -f2`",
				"Rscript batch-consistency-analysis.r \$FILE1 \$FILE2 -1 \$OUTPUT 0 F signal.value"
			);
		}
		if ($Config{'peak'} eq "macs2" || $Config{'peak'} eq "both"){
			printf OUT2 ("\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n",
				"#macs2",
				"FILE3=`echo FILES | cut -d\" \" -f1 | sed \'s\/\.regionPeak\.gz\/\.narrowPeak\.gz\/\'`",
				"FILE4=`echo FILES | cut -d\" \" -f2 | sed \'s\/\.regionPeak\.gz\/\.narrowPeak\.gz\/\'`",
				"Rscript batch-consistency-analysis.r \$FILE3 \$FILE4 -1 \$\{OUTPUT\}.macs 0 F p.value"
			);
		}
		printf OUT2 ("\t\t%s\n\t%s\n\t\t%s\n\t\t%s\n",
			"\'\;",
			"# Generate bed files",
			"max_numPeaks_Rep_spp=0",
			"max_numPeaks_Rep_macs=0"
		);
                for (my$i=0;$i<@IPs;$i++){
                        my $R1=substr($IPs[$i],-1);
                        print OUT2 "\t\t## IDR with Rep$R1\n";
			if ($Config{'peak'} eq "spp" || $Config{'peak'} eq "both"){
				print OUT2 "\t\t\tnumPeaks_Rep$R1\_pr_spp=\$(awk '\$11 <= $IDR_pr {print \$0}' $IPs[$i].pr1_VS_$IPs[$i].pr2-overlapped-peaks.txt | wc -l)\n";
			}
			if ($Config{'peak'} eq "macs2" || $Config{'peak'} eq "both"){
				print OUT2 "\t\t\tnumPeaks_Rep$R1\_pr_macs=\$(awk '\$11 <= $IDR_pr {print \$0}' $IPs[$i].pr1_VS_$IPs[$i].pr2.macs-overlapped-peaks.txt | wc -l)\n";
			}
			for (my$j=$i+1;$j<@IPs;$j++){
                                my$R2=substr($IPs[$j],-1);
				if ($Config{'peak'} eq "spp" || $Config{'peak'} eq "both"){
					print OUT2 "\t\t\tnumPeaks_Rep$R2\_Rep$R1\_spp=\$(awk '\$11 <= $IDR_true {print \$0}' $IPs[$j]_VS_$IPs[$i]-overlapped-peaks.txt | wc -l)\n";
					print OUT2 "\t\t\tif [ \$numPeaks_Rep$R2\_Rep$R1\_spp -gt \$max_numPeaks_Rep_spp ];\n\t\t\t\tthen max_numPeaks_Rep_spp=\$numPeaks_Rep$R2\_Rep$R1\_spp\n\t\t\tfi\n";
				}
				if ($Config{'peak'} eq "macs2" || $Config{'peak'} eq "both"){
					print OUT2 "\t\t\tnumPeaks_Rep$R2\_Rep$R1\_macs=\$(awk '\$11 <= $IDR_true {print \$0}' $IPs[$j]_VS_$IPs[$i].macs-overlapped-peaks.txt | wc -l)\n";
					print OUT2 "\t\t\tif [ \$numPeaks_Rep$R2\_Rep$R1\_macs -gt \$max_numPeaks_Rep_macs ];\n\t\t\t\tthen max_numPeaks_Rep_macs=\$numPeaks_Rep$R2\_Rep$R1\_macs\n\t\t\tfi\n";
				}
                        }
                }
                print OUT2 "\t\t## IDR for combined Rep\n";
		if ($Config{'peak'} eq "spp" || $Config{'peak'} eq "both"){
			print OUT2 "\t\t\tnumPeaks_Rep0_pr_spp=\$(awk '\$11 <= $IDR_pr {print \$0}' $d\_IP_Rep0.pr1_VS_$d\_IP_Rep0.pr2-overlapped-peaks.txt | wc -l)\n";
		}
		if ($Config{'peak'} eq "macs2" || $Config{'peak'} eq "both"){
			print OUT2 "\t\t\tnumPeaks_Rep0_pr_macs=\$(awk '\$11 <= $IDR_pr {print \$0}' $d\_IP_Rep0.pr1_VS_$d\_IP_Rep0.pr2.macs-overlapped-peaks.txt | wc -l)\n";
		}
                print OUT2 "\t\t## print optimal and conservitive IDR values\n";
		if ($Config{'peak'} eq "spp" || $Config{'peak'} eq "both"){
			printf OUT2 ("\t\t\t%s\n\t\t\t\t%s\n\t\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n",
				"if [ \$numPeaks_Rep0_pr_spp -gt \$max_numPeaks_Rep_spp ];",
				"then optThresh_spp=\$numPeaks_Rep0_pr_spp",
				"else optThresh_spp=\$max_numPeaks_Rep_spp",
				"fi",
				"zcat $d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz | sort -k7nr,7nr | head -n \$optThresh_spp | pigz -c > spp.optimal.$d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz",
				"zcat $d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz | sort -k7nr,7nr | head -n \$max_numPeaks_Rep_spp | pigz -c >  spp.conservative.$d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.regionPeak.gz"
			);
		}
		if ($Config{'peak'} eq "macs2" || $Config{'peak'} eq "both"){
			printf OUT2 ("\t\t\t%s\n\t\t\t\t%s\n\t\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n",
				"if [ \$numPeaks_Rep0_pr_macs -gt \$max_numPeaks_Rep_macs ];",
				"then optThresh_macs=\$numPeaks_Rep0_pr_macs",
				"else optThresh_macs=\$max_numPeaks_Rep_macs",
				"fi",
				"zcat $d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.narrowPeak.gz | sort -k7nr,7nr | head -n \$optThresh\_macs | pigz -c > macs.optimal.$d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.narrowPeak.gz",
				"zcat $d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.narrowPeak.gz | sort -k7nr,7nr | head -n \$max_numPeaks_Rep\_macs | pigz -c >  macs.conservative.$d\_IP_Rep0.tagAlign_VS_$d\_Input_Rep0.tagAlign.narrowPeak.gz"
			);
		}
		# print lines for report.sh
		printf OUT3 ("\t%s\n\t%s\n",
			"max_numPeaks_Rep_spp=0",
			"max_numPeaks_Rep_macs=0"
		);
		for (my$i=0;$i<@IPs;$i++){
			my $R1=substr($IPs[$i],-1);
			printf OUT3 ("\t%s\n\t\t%s\n\t\t%s\n\t\t%s\n\t\t\t%s\n\t\t%s\n\t\t%s\n\t\t\t%s\n\t\t%s\n",
				"## IDR with Rep$R1",
				"numPeaks_Rep$R1\_pr_spp=0",
				"numPeaks_Rep$R1\_pr_macs=0",
				"if \[ -e $IPs[$i].pr1_VS_$IPs[$i].pr2-overlapped-peaks.txt \]; then",
					"numPeaks_Rep$R1\_pr_spp=\$(awk '\$11 <= $IDR_pr {print \$0}' $IPs[$i].pr1_VS_$IPs[$i].pr2-overlapped-peaks.txt | wc -l)",
				"fi",
				"if \[ -e $IPs[$i].pr1_VS_$IPs[$i].pr2.macs-overlapped-peaks.txt \]; then",
					"numPeaks_Rep$R1\_pr_macs=\$(awk '\$11 <= $IDR_pr {print \$0}' $IPs[$i].pr1_VS_$IPs[$i].pr2.macs-overlapped-peaks.txt | wc -l)",
				"fi"
			);
			if ($Config{'peak'} eq "spp"){
				print OUT3 "\t\techo \"numPeaks_Rep$R1\_pr\t\$numPeaks_Rep$R1\_pr_spp\tN\/A\"\n";
			}
			if ($Config{'peak'} eq "macs2"){
				print OUT3 "\t\techo \"numPeaks_Rep$R1\_pr\tN\/A\t\$numPeaks_Rep$R1\_pr_macs\"\n";
			}
			if ($Config{'peak'} eq "both"){
				print OUT3 "\t\techo \"numPeaks_Rep$R1\_pr\t\$numPeaks_Rep$R1\_pr_spp\t\$numPeaks_Rep$R1\_pr_macs\"\n";
			}
			for (my$j=$i+1;$j<@IPs;$j++){
				my$R2=substr($IPs[$j],-1);
				printf OUT3 ("\n\t\t%s\n\t\t%s\n\t\t%s\n\t\t\t%s\n\t\t%s\n\t\t%s\n\t\t\t%s\n\t\t%s\n\t\t%s\n\t\t\t%s\n\t\t%s\n\t\t%s\n\t\t\t%s\n\t\t%s\n",
					"numPeaks_Rep$R2\_Rep$R1\_spp=0",
					"numPeaks_Rep$R2\_Rep$R1\_macs=0",
					"if \[ -e $IPs[$j]_VS_$IPs[$i]-overlapped-peaks.txt \]; then",
						"numPeaks_Rep$R2\_Rep$R1\_spp=\$(awk '\$11 <= $IDR_true {print \$0}' $IPs[$j]_VS_$IPs[$i]-overlapped-peaks.txt | wc -l)",
					"fi",
					"if \[ -e $IPs[$j]_VS_$IPs[$i].macs-overlapped-peaks.txt \]; then",
						"numPeaks_Rep$R2\_Rep$R1\_macs=\$(awk '\$11 <= $IDR_true {print \$0}' $IPs[$j]_VS_$IPs[$i]\.macs-overlapped-peaks.txt | wc -l)",
					"fi",
					"if [ \$numPeaks_Rep$R2\_Rep$R1\_spp -gt \$max_numPeaks_Rep_spp ];",
						"then max_numPeaks_Rep_spp=\$numPeaks_Rep$R2\_Rep$R1\_spp",
					"fi",
					"if [ \$numPeaks_Rep$R2\_Rep$R1\_macs -gt \$max_numPeaks_Rep_macs ];",
						"then max_numPeaks_Rep_macs=\$numPeaks_Rep$R2\_Rep$R1\_macs",
					"fi"
				);
				if ($Config{'peak'} eq "spp"){
					print OUT3 "\t\techo \"numPeaks_Rep$R2\_Rep$R1\t\$numPeaks_Rep$R2\_Rep$R1\_spp\tN\/A\"\n";
				}
				if ($Config{'peak'} eq "macs2"){
					print OUT3 "\t\techo \"numPeaks_Rep$R2\_Rep$R1\tN\/A\t\$numPeaks_Rep$R2\_Rep$R1\_macs\"\n";
				}
				if ($Config{'peak'} eq "both"){
					print OUT3 "\t\techo \"numPeaks_Rep$R2\_Rep$R1\t\$numPeaks_Rep$R2\_Rep$R1\_spp\t\$numPeaks_Rep$R2\_Rep$R1\_macs\"\n";
				}
			}
		}
		printf OUT3 ("\t%s\n\t\t%s\n\t\t%s\n\t\t%s\n\t\t\t%s\n\t\t%s\n\t\t%s\n\t\t\t%s\n\t\t%s\n",
			"## IDR for combined Rep",
			"numPeaks_Rep0_pr_spp=0",
			"numPeaks_Rep0_pr_macs=0",
			"if \[ -e $d\_IP_Rep0.pr1_VS_$d\_IP_Rep0.pr2-overlapped-peaks.txt \]; then",
				"numPeaks_Rep0_pr_spp=\$(awk '\$11 <= $IDR_pr {print \$0}' $d\_IP_Rep0.pr1_VS_$d\_IP_Rep0.pr2-overlapped-peaks.txt | wc -l)",
			"fi",
			"if \[ -e $d\_IP_Rep0.pr1_VS_$d\_IP_Rep0.pr2.macs-overlapped-peaks.txt \]; then",
				"numPeaks_Rep0_pr_macs=\$(awk '\$11 <= $IDR_pr {print \$0}' $d\_IP_Rep0.pr1_VS_$d\_IP_Rep0.pr2.macs-overlapped-peaks.txt | wc -l)",
			"fi",
		);
		if ($Config{'peak'} eq "spp"){
			print OUT3 "\t\techo \"numPeaks_Rep0_pr\t\$numPeaks\_Rep0\_pr_spp\tN\/A\"\n";
		}
		if ($Config{'peak'} eq "macs2"){
			print OUT3 "\t\techo \"numPeaks_Rep0_pr\tN\/A\t\$numPeaks\_Rep0\_pr\_macs\"\n";
		}
		if ($Config{'peak'} eq "both"){
			print OUT3 "\t\techo \"numPeaks_Rep0_pr\t\$numPeaks\_Rep0\_pr_spp\t\$numPeaks\_Rep0\_pr\_macs\"\n";
		}
		printf OUT3 ("\t%s\n\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t%s\n\t\t%s\n\t\t\t%s\n\t\t\t%s\n\t\t%s\n",
			"## print optimal and conservitive IDR values",
			"if [ \$numPeaks_Rep0_pr_spp -gt \$max_numPeaks_Rep_spp ];",
				"then optThresh_spp=\$numPeaks_Rep0_pr_spp",
				"else optThresh_spp=\$max_numPeaks_Rep_spp",
			"fi",
			"if [ \$numPeaks_Rep0_pr_macs -gt \$max_numPeaks_Rep_macs ];",
				"then optThresh_macs=\$numPeaks_Rep0_pr_macs",
				"else optThresh_macs=\$max_numPeaks_Rep_macs",
			"fi"
		);
		if ($Config{'peak'} eq "spp"){
			printf OUT3 ("\t\t%s\n\t\t%s",
				"echo \"optThresh\t\$optThresh_spp\tN\/A\"",
				"echo \"conThresh\t\$max_numPeaks_Rep_spp\tN\/A\"\n"
			);
		}
		if ($Config{'peak'} eq "macs2"){
			printf OUT3 ("\t\t%s\n\t\t%s",
				"echo \"optThresh\tN\/A\t\$optThresh_macs\"",
				"echo \"conThresh\tN\/A\t\$max_numPeaks_Rep_macs\"\n"
			);
		}
		if ($Config{'peak'} eq "both"){
			printf OUT3 ("\t\t%s\n\t\t%s",
				"echo \"optThresh\t\$optThresh_spp\t\$optThresh_macs\"",
				"echo \"conThresh\t\$max_numPeaks_Rep_spp\t\$max_numPeaks_Rep_macs\"\n"
			);
		}
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
ls *narrowPeak.gz >> regionPeak.list
while read -r line; do
	case \$line in \*\"regionPeak.gz\"\)
		zcat \$line | cut -f1-3 > \$line.bed33
	esac
	case \$line in \*\"narrowPeak.gz\"\)
		zcat \$line | cut -f1-3 > \$line.bed33
	esac
done < regionPeak.list
mv *narrowPeak.gz.bed33 bed/
mv *regionPeak.gz.bed33 bed/
rename -f s/regionPeak.gz.bed33/bed/ bed/*
rename -f s/narrowPeak.gz.bed33/macs.bed/ bed/*

### removed unzipped files
rm -f *tagAlign

### remove 50% overlapping peaks from blacklist
ls spp.*.regionPeak.gz  > spp.list
while read -r line; do
	bedtools intersect -a \$line -b $bin_dir/$black_list  -f 0.5 -v > \$line.rmblacklist
done < spp.list
rename s/regionPeak.gz.rmblacklist/regionPeak.rmblacklist/ *.regionPeak.gz.rmblacklist
ls macs.*.narrowPeak.gz > macs.list
while read -r line; do
	bedtools subtract -a \$line -b $bin_dir/$black_list > \$line.rmblacklist
done < macs.list
rename s/narrowPeak.gz.rmblacklist/narrowPeak.rmblacklist/ *.narrowPeak.gz.rmblacklist

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
