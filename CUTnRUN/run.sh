for file in SRR*.gz
#for file in SRR5509755.fastq.gz
do
  echo "running "$file
  pigz -dvf $file;
  wait;
  ./STAR/bin/Linux_x86_64/STAR --genomeDir ./STAR_GENOME --runThreadN 8 --outSAMstrandField intronMotif --sjdbGTFfile c_elegans.PRJNA13758.WS267.canonical_geneset.gtf --readFilesIn ${file%.gz} --outFileNamePrefix ${file%fastq.gz} --outSAMtype BAM SortedByCoordinate
  cufflinks -p 16 -g c_elegans.PRJNA13758.WS267.canonical_geneset.gtf ${file%fastq.gz}Aligned.sortedByCoord.out.bam;
  mv isoforms.fpkm_tracking ${file%fastq.gz}fpmk_tracking;
  mv transcripts.gtf ${file%fastq.gz}transcripts.gtf;
  cat ${file%.gz} | pigz -c > $file;
  rm ${file%.gz};
done
