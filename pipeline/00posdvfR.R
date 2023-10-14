setwd("/Users/02dvf/")

df<-read.delim("../sample-name-all.txt",header=FALSE,check.names = FALSE) 
df<-as.data.frame(df)
library(magrittr) 
library(dplyr)

for(i in 1:nrow(df))
{
  #i=1
  oldname = paste("./",df[i,],"_contigs.fasta_gt10000bp_dvfpred.txt",sep='')
  result = read.delim(oldname,sep='\t', header=FALSE,dec = ".")
  names(result)<-c("name","len","score","pvalue")
  new.data <- result[which(result$score>0.9 & result$pvalue<0.05),] #去除len小于1000
  newname = paste("../04posdvf/",df[i,],".contigs0.9-0.05-10kdvf.txt",sep='')
  write.table(new.data,newname,row.names = FALSE,col.names = FALSE,quote =FALSE,sep='\t')
  new.data_name <- as.data.frame(new.data[,1])
  newname_name = paste("../04posdvf/",df[i,],".contigs0.9-0.05-10kdvf_name.txt",sep='')
  write.table(new.data_name,newname_name,row.names = FALSE,col.names = FALSE,quote =FALSE,sep='\t')
}




