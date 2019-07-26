files <- list.files(path="./",pattern="*.lengths$",recursive=FALSE)
library(ggplot2)

for (i in 1:length(files)){
  print(files[i])
  rdata <- read.delim(files[i], header=TRUE,sep="\t")
  p=ggplot(data=rdata, aes(x=rdata[,1],y=rdata[,2])) + geom_bar(stat="identity", fill="black",width=1) + labs(x="Insert Lengths (bp)",y="Count",title=paste("Length Distrobution",files[i]))
  ggsave(p,file=paste(files[i],".png",sep=""))
}
