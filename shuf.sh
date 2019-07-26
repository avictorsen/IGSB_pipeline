#!/bin/sh
fileName=$1
outputDir=$2
outputStub=$3
#fileName='Dmel_WPP_CTCF_pooled.txt.map.sort.bam.bed' # input tagAlign file name
#outputDir='/glusterfs/users/malijia/ComparativeTF/data/CTCF/WPP/PooledPrep/' # output directory for pseudoreplicate files
#outputStub='Dmel_WPP_CTCF_pooled' # prefix name for pseudoReplicate files

echo ""
echo "file: $1"
echo "outdir: $2"
echo "outprefix: $3"
echo "zcat ${fileName} | wc -l ";
nlines=$( zcat ${fileName} | wc -l ) # Number of reads in the tagAlign file
#echo "((nlines + 1) / 2)";
nlines=$(( (nlines + 1) / 2 )) # half that number
echo "zcat \"${fileName}\" | shuf | split -d -l ${nlines} - \"${outputDir}/${outputStub}\"";
zcat "${fileName}" | shuf | split -d -l ${nlines} - "${outputDir}/${outputStub}" # This will shuffle the lines in the file and split it into two parts
#gzip "${outputDir}/${outputStub}00"
#gzip "${outputDir}/${outputStub}01"
pigz "${outputDir}/${outputStub}00"
pigz "${outputDir}/${outputStub}01"
echo "mv \"${outputDir}/${outputStub}00.gz\" \"${outputDir}/${outputStub}.pr1.tagAlign.gz\"";
mv "${outputDir}/${outputStub}00.gz" "${outputDir}/${outputStub}.pr1.tagAlign.gz"
echo "mv \"${outputDir}/${outputStub}01.gz\" \"${outputDir}/${outputStub}.pr2.tagAlign.gz\"";
mv "${outputDir}/${outputStub}01.gz" "${outputDir}/${outputStub}.pr2.tagAlign.gz"
wait
