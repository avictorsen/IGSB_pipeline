rm(list=ls())

# ##
# Parantu Shah, April 2009
# ##

# ##
# Specify the directory, input and control file.
# ##

args <- commandArgs()

baseDir <- sub('--baseDir=', '', args[grep('--baseDir=', args)])
outDir <- sub('--outDir=', '', args[grep('--outDir=', args)])
factor <- sub('--factor=', '', args[grep('--factor=', args)])
ipFileName <- sub('--ipFile=', '', args[grep('--ipFile=', args)])
controlFileName <- sub('--controlFile=', '', args[grep('--controlFile=', args)])


ipFile <- paste (baseDir, ipFileName, sep = "/")
controlFile <- paste (baseDir, controlFileName, sep = "/")

pdfString <- paste (factor, "_crosscorrelation.pdf", sep="")
pdfFile <- paste (outDir, pdfString, sep="/")

densityWigString <- paste (factor, "_density.wig", sep="")
densityFile <- paste (outDir, densityWigString, sep="/")   		
wigString <- paste (factor, ", smoothed, background-subtracted tag density", sep=" ")

enrichmentWigString <- paste (factor, "_enrichment.estimates.wig", sep="")
enrichmentFile <- paste (outDir, enrichmentWigString, sep="/")
wigEnricmentString <- paste (factor, ", conservative fold-enrichment/depletion estimates shown on log2 scale", sep=" ")

# ##

library(spp)
library(snow)

#cluster <- makeCluster(8);

chip.data <- read.bam.tags(ipFile);
input.data <- read.bam.tags(controlFile);

binding.characteristics <- get.binding.characteristics(chip.data,srange=c(50,500),bin=5);
#print(paste("binding peak separation distance =",binding.characteristics$peak$x))

pdf(file=pdfFile,width=5,height=5)
par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8);
plot(binding.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation");
abline(v=binding.characteristics$peak$x,lty=2,col=2)
dev.off();

chip.data <- select.informative.tags(chip.data,binding.characteristics);
input.data <- select.informative.tags(input.data,binding.characteristics);

#chip.data <- remove.local.tag.anomalies(chip.data);
#input.data <- remove.local.tag.anomalies(input.data);

tag.shift <- round(binding.characteristics$peak$x/2)
smoothed.density <- get.smoothed.tag.density(chip.data,control.tags=input.data,bandwidth=200,step=100,tag.shift=tag.shift);

writewig(smoothed.density, densityFile, wigString);
rm(smoothed.density);

enrichment.estimates <- get.conservative.fold.enrichment.profile(chip.data,input.data,fws=500,step=100,alpha=0.01);
writewig(enrichment.estimates, enrichmentFile, wigEnricmentString);
rm(enrichment.estimates);

