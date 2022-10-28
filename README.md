Project name : Gene  expression analysis and survival analysis

Introduction:The determination of the pattern of genes expressed at the level of genetic transcription, under specific circumstances or in a specific cell,is called as gene expression analysis.The expression studies are directed to detect and quantify messenger RNA (mRNA) levels of a specific gene.It  comprises from the gene activation until the mature protein is located in its corresponding compartment to perform its function and contribute to the expression of the phenotype of cell.

 Here I choose Cervical cancer dataset downloaded from NCBI GEO(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168009)

1.Normalization:It is a way of organizing a dataset .It removes the redundancy and unstructured data from the dataset.
Thus the dataset was normalized as given below:

cpm_matrix[,i]=(data_file[,i]/sum(data_file[,i]))*1000000

logcpm=log2(cpm_matrix+1)

2.Heatmap:It visualizes the graphical representation of the data using colors,here in the analysis part the z-scores were being visualized in the mentioned dataset.

3.Differential Gene analysis:It takes the read count data and performs the statistical analysis to discover the quantivative changes between experimental groups.
Here I categorized the data assuming NDB1,NDB2,NDB3,NDB4 as control samples and  DCB1,DCB2,DCB3,DCB4,DCB5 as tumor samples.
log2fc of the normalized data was calculated,it reflects how different the expression of a gene in one condition is from the expression of the same gene in another condition.

vec1=gene[i,control]

vec2=gene[i,tumor]

log2cpm=mean(vec1) - mean(vec2)

4.Circos plot: It represents information with long axes or a large amount of categories; second, it intuitively shows data with multiple tracks focusing on the same object; third, it easily demonstrates relations between elements and the plots are in the circular format.

5.SSGSEA:Stands for "Single Sample Geneset Enrichment Analysis ".It calculates separate enrichment scores for each pairing of a sample and gene set. Each ssGSEA enrichment score represents the degree to which the genes in a particular gene set are coordinately up- or down-regulated within a sample.

6.Survival analysis:It involves statistical measures and implementations for data analysis where the outcome variable of interest is time until an event occurs.

7.Single cell RNA-Seq analysis:scRNA allows us to understand cellular differences in expression, and hence it is directly applicable to the studies of cell heterogeneity, cell population and subpopulation identification, effects of low copy mRNA distribution and transcriptional regulation.Single cell RNA-seq data used here is generated using the 10X genomics platform.

Here I took cervical cancer dataset for single cell RNA-Seq analysis(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190075) and from the mentioned dataset I choose the files given below

GSM5712850_CC04_barcodes.tsv.gz

GSM5712850_CC04_features.tsv.gz

GSM5712850_CC04_matrix.mtx.gz

as custom files for Single cell RNA-Seq analysis,renamed them as barcodes.tsv.gz,features.tsv.gz,matrix.mtx.gz .
It further visualizes the violin plot,elbow plot,heatmaps,scatterplots.

Name:Khushbu Narendra Sharma
