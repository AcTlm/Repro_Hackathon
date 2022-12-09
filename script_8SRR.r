## AUTHOR: AC TOULEMONDE
## DATE = 20/11/2022

#######  PACKAGES 

#install.packages("RColorBrewer")
#install.packages("ggplot2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library( "DESeq2")
library("ggplot2")
library("readr")
library("RColorBrewer")


#### IMPORT DATA AND TRANSFORM
setwd("~/AMI2B/Hackathon/DESeq2_analyse_stat/Rendu_hackathon/plot")

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
a = plotDispEsts(dds, main = "Dispersion des gènes")
dev.copy2pdf(file = "dispersion_genes.pdf")

#### Analyse of replication by group 
rld <- rlog(dds)  ## log of normalized counts 
PCA = plotPCA(rld,intgroup = "cellType") ## Principal Components Analyise by group cancer & human_ref ? 
PCA + geom_label(aes(label = sampleName))
PCA + ggtitle ("Principal Components analyis plot")
dev.copy2pdf(file = "PCA_Hackathon.pdf")


#### Heatmap 
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
dists <- dist(t(assay(rld)))
mat <- as.matrix(dists)
rownames(mat) <- colnames(mat) <- as.character(colData$sampleName)
hc <- hclust(dists)

heatmap(mat, Rowv=as.dendrogram(hc), symm=T, trace='none', col=rev(hmcol), margin=c(13,13), main='Heatmap')
dev.copy2pdf(file = "heatmap.pdf")

#### differential expression between treatment conditions 
res <- results(dds)#,contrast = c("cellType","cancer","human"),alpha = 0.05)
result = res[order(res$padj),]
result = as.data.frame(result)
head(result)

hist( res$pvalue, breaks=20, col="grey" )
dev.copy2pdf(file = "residuals_pvalue.pdf")


#### log2Fold Changes as function of expression----------------
plotMA(results(dds),alpha = 0.05, )
abline(h=c(-2,2),col = "blue")

#### Volcano plot sur les genes etudiées en biblio 
result_analyse<-result[c("ENSG00000088256","ENSG00000101019","ENSG00000114770","ENSG00000115524","ENSG00000156052","ENSG00000163930","ENSG00000245694","ENSG00000148848"),]
head(result_analyse)

par(mfrow=c(1,2))
a = -log10(0.05)

# Make volcano plot 
### logfold2Change exprié 2/4 fois plus que les autres 
with(result_analyse, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot on studied genes", xlim=c(-10,10)))
with(subset(result_analyse, pvalue<.05 & log2FoldChange>2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(result_analyse, pvalue<.05 & log2FoldChange<(-2)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
abline(h = a) 
abline(v = (-2))
abline(v = (2))
dev.copy2pdf(file = "volcano_plot_genes_studied.pdf")

# Make a basic volcano plot 
with(result, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot for all genes ", xlim=c(-10,10)))
with(subset(result, pvalue<.05 & log2FoldChange>2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(result, pvalue<.05 & log2FoldChange<(-2)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
abline(h = a)
abline(v = (-2))
abline(v = (2))
dev.copy2pdf(file = "volcano_plot_all_genes.pdf")

##### 
#setwd("~/AMI2B/Hackathon/DESeq2_analyse_stat/Rendu_hackathon")
#write.table(result_p,file="resultats_DESEQ_pvalue.csv",sep =" ")


