rm(list=ls())
args <- commandArgs()
print(args) 
baseDir <- sub('--baseDir=', '', args[grep('--baseDir=', args)])
outDir <- sub('--outDir=', '', args[grep('--outDir=', args)])
factor <- sub('--factor=', '', args[grep('--factor=', args)])
ipFileName <- sub('--ipFile=', '', args[grep('--ipFile=', args)])

ipFile <- paste (baseDir, ipFileName, sep = "/")

densityWigString1 <- paste (factor, "_density.wig", sep="")
densityWigString2 <- paste (factor, "_density.scaled.wig", sep="")
densityFile1 <- paste (outDir, densityWigString1, sep="/")
densityFile2 <- paste (outDir, densityWigString2, sep="/")
wigString1 <- paste (factor, ", smoothed tag density", sep="")
wigString2 <- paste (factor, ", smoothed tag density, scaled by dataset size", sep="")

# ##

library(spp)

#chip.data <- read.tagalign.tags(ipFile);
chip.data <- read.bam.tags(ipFile);
binding.characteristics <- get.binding.characteristics(chip.data,srange=c(50,500),bin=5);

chip.data <- select.informative.tags(chip.data,binding.characteristics);

#tag.shift <- round(binding.characteristics$peak$x/2)

# Modified to scale the results by dataset size

#smoothed.density <- get.smoothed.tag.density(chip.data,bandwidth=200,step=100,tag.shift=tag.shift);
smoothed.density <- get.smoothed.tag.density(chip.data,bandwidth=200,step=100);
writewig(smoothed.density, densityFile1, wigString1);
rm(smoothed.density);

#smoothed.density <- get.smoothed.tag.density(chip.data,bandwidth=200,step=100,tag.shift=tag.shift,scale.by.dataset.size=T);
smoothed.density <- get.smoothed.tag.density(chip.data,bandwidth=200,step=100,scale.by.dataset.size=T);
writewig(smoothed.density, densityFile2, wigString2);
rm(smoothed.density);

