setwd("/Users/02visorter/")
df<-read.delim("../sample-name-all.txt",header=FALSE,check.names = FALSE) 
df<-as.data.frame(df)
library(magrittr) 
library(dplyr)

for(i in 1:nrow(df))
{
  oldname = paste("./",df[i,],"final-viral-score.tsv",sep='')
  new.data = read.delim(oldname,header=TRUE,dec = ".")
  new.data <- as.data.frame(new.data$seqname)
  names(new.data)<-c("seqname")
  newname = paste("../06vs2name/",df[i,],"viralseqname.txt",sep='')
  write.table(new.data,newname,row.names = FALSE,col.names = TRUE,quote =FALSE)
}

