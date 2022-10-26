#Cancer dataset
setwd('C:\\Users\\sai\\Documents\\R')
data_file=read.csv("GSE168009_Raw_count.csv",row.names=1)
#count per matrix
cpm_matrix=data_file
for(i in 1:ncol(data_file)){
  cpm_matrix[,i]=(data_file[,i]/sum(data_file[,i]))*1000000
}
cpm_matrix[is.na(cpm_matrix)]=0

#log of cpm
logcpm=log2(cpm_matrix+1)
summary(logcpm)



