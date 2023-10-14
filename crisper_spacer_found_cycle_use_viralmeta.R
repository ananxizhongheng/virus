setwd("/users/spades/")
crispr_file <- dir('01crisper')
for (file in crispr_file) {
  data <- read.delim(paste('01crisper', file, sep = '/'), sep = '\t', header=F,stringsAsFactors = FALSE, check.names = FALSE)
  
  data_spacer=as.data.frame(data[-c(1:6),4])
  names(data_spacer)<-c("node_name")
  
  data_spacer<-as.data.frame(data_spacer[!grepl("Average", data_spacer$node_name),]) #去除含有average的行
  
  names(data_spacer)<-c("node_name")
  
  df <- as.data.frame(data_spacer[order(data_spacer$node_name,decreasing = T),]) #降序排列
  names(df)<-c("node_name")
  df<-as.data.frame(df[!apply(df == "", 1, all),]) #去除空白行
  write.table(df,paste("./01crisper/",file,"viralcontigs-crt-clean.txt",sep = ''),row.names = F,col.names =F, quote =F)
}



