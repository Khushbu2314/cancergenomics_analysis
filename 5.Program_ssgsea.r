#cancer dataset
setwd('C:\\Users\\sai\\Documents\\R')
data_file=read.csv("GSE168009_Raw_count.csv",row.names = 1)

#Create a count per matrix
cpm_matrix=data_file
for(i in 1:ncol(data_file)){
  cpm_matrix[,i]=(data_file[,i]/sum(data_file[,i]))*1000000
}
cpm_matrix[is.na(cpm_matrix)]=0

#log of cpm
logcpm=log2(cpm_matrix+1)
summary(logcpm)

library(matrixStats)
library(circlize)
library(ComplexHeatmap)
#library(data.table)

ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

saveRDS(logcpm,file="input1.rds")
data = readRDS("input1.rds")
data[1:5,1:5]
data = as.matrix(data)

gene_set = read.csv("markers2Sep.csv")
head(gene_set)

gene_sets = as.list(as.data.frame(gene_set))
print("genes set ready")

res=ssgsea(data,gene_set,norm = F,scale = T)
print(res)

