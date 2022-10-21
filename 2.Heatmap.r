data_file=read.csv("GSE168009_Raw_count.csv",row.names=1)
#z-score
library(matrixStats)
z_score = (logcpm - rowMeans(logcpm))/rowSds(as.matrix(logcpm))[row(logcpm)]
z_score[is.na(z_score)]=0
print(z_score)

#Variance
variance = apply(logcpm, 1, var)
variance = sort(variance,decreasing = T)
top50 = variance[1:50]
pmat = z_score[names(top50),]


#heatmap
library(ComplexHeatmap)
Heatmap(pmat)
