## AUTHOR: AC TOULEMONDE
## DATE = 20/11/2022

#######  PACKAGES 

install.packages("RColorBrewer")
install.packages("ggfortify")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library( "DESeq2")
library("ggplot2")
library("readr")
library("RColorBrewer")
library("ggfortify")

#### IMPORT DATA AND TRANSFORM
setwd("~/AMI2B/Hackathon/DESeq2_analyse_stat")
#counts_db= read_delim("count_8SRR.csv", delim = ";", escape_double = FALSE, col_types = cols(SRR628582.bam = col_integer(),SRR628583.bam = col_integer(), SRR628584.bam = col_integer(),SRR628585.bam = col_integer(), SRR628586.bam = col_integer(),SRR628587.bam = col_integer(), SRR628588.bam = col_integer(), SRR628589.bam = col_integer()), trim_ws = TRUE)

counts_db <- read.delim("~/AMI2B/Hackathon/DESeq2_analyse_stat/Rendu_hackathon/count.txt", comment.char="#", stringsAsFactors=TRUE)
annotation <- counts_db[,c(1:6)]
head(annotation)
rownames(counts_db) <- counts_db$Geneid
counts <- counts_db[,-c(1:6)]
colnames(counts) <- gsub("Alignments\\.|_sorted.bam","",colnames(counts))


cellType <- c(rep("Mutant",3),rep("wild",5))
sampleName <- paste(cellType,c(1:3),sep = "_")
colData <- data.frame(row.names = colnames(counts), cellType = cellType,sampleName = sampleName)
colData[,1] = as.factor(colData[,1])

dds <- DESeqDataSetFromMatrix(countData = counts, design = ~ cellType, colData = colData)
mcols(dds) <- DataFrame(annotation)

#### normalizes the read counts,estimates dispersions, and fits the linear model
dds <- DESeq(dds,betaPrior = T)
sizeFactors(dds) ## normalization factors


####  Dispersion plot
plotDispEsts(dds)
#dev.copy2pdf(file = "dispEsts.pdf")

#### Analyse of replication by group 
rld <- rlog(dds)  ## log of normalized counts 
PCA = plotPCA(rld,intgroup = "cellType") ## Principal Components Analyise by group cancer & human_ref ? 
PCA + geom_label(aes(label = sampleName))
#dev.copy2pdf(file = "PCA_Hackathon")


#### Heatmap 
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
dists <- dist(t(assay(rld)))
mat <- as.matrix(dists)
rownames(mat) <- colnames(mat) <- as.character(colData$sampleName)
hc <- hclust(dists)

heatmap(mat, Rowv=as.dendrogram(hc), symm=T, trace='none', col=rev(hmcol), margin=c(13,13), main='heatmap')
#dev.copy2pdf(file = "distClustering.pdf")

#### differential expression between treatment conditions 
res <- results(dds)#,contrast = c("cellType","cancer","human"),alpha = 0.05)
result_p = res[order(res$padj),]
result_p = as.data.frame(result_p)

summary(res)
#setwd("~/AMI2B/Hackathon/DESeq2_analyse_stat/Rendu_hackathon")
#write.table(result_p,file="resultats_DESEQ_pvalue.csv",sep =" ")

hist( res$pvalue, breaks=20, col="grey" )

a = plotCounts(dds,gene = "Id4",intgroup = "cellType")
dev.copy2pdf(file = "Id4_plotCounts.pdf")
subset(res,grepl("Hox",rownames(res))) ## subset the results table to look at your favourite genes:

#### log2 Fold Changes as function of expression----------------
plotMA(results(dds),alpha = 0.05)
abline(h=c(-2,2),col = "blue")


setwd("~/AMI2B/Hackathon/DESeq2_analyse_stat/Rendu_hackathon")
write.table(result_p,file="resultats_DESEQ_pvalue.csv",sep =" ")


