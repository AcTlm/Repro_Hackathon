
library( "DESeq2")
library("ggplot2")
library("readr")
library("RColorBrewer")
library("fdrtool")


### Importation et structiration des données 
counts_db <- read.delim("~/AMI2B/Hackathon/DESeq2_analyse_stat/Rendu_hackathon/count.txt", comment.char="#", stringsAsFactors=TRUE)
counts_db <- read.delim("~/AMI2B/Hackathon/DESeq2_analyse_stat/Rendu_hackathon/count.txt", comment.char="#", stringsAsFactors=TRUE)
annotation <- counts_db[,c(1:6)]
rownames(counts_db) <- counts_db$Geneid
counts <- counts_db[,-c(1:6)]
colnames(counts) <- gsub("Alignments\\.|_sorted.bam","",colnames(counts))


cellType <- c(rep("Mutant",3),rep("wild",5))
sampleName <- paste(cellType,c(1:3),sep = "_")
colData <- data.frame(row.names = colnames(counts), cellType = cellType,sampleName = sampleName)
colData[,1] = as.factor(colData[,1])

### Construction des métadonnées
dds <- DESeqDataSetFromMatrix(countData = counts, design = ~ cellType, colData = colData)
mcols(dds) <- DataFrame(annotation)


###  Conrole qualité et normalisation des données de comptage
GeneCounts <- counts(dds)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)
nz.counts <- subset(GeneCounts, idx.nz)
sam <- sample(dim(nz.counts)[1], 5)
nz.counts[sam, ]


###  Normalisation et estimation des facteurs de taill des SRR
dds <- DESeq(dds,betaPrior = T)
sizeFactors(dds) ## normalization factors
plotDispEsts(dds, main = "Dispersion des gènes")
dev.copy2pdf(file = "dispersion_genes.pdf")

####PCA et cartes thermiques
rld = rlogTransformation(dds)
##  Heatmap 
dists <- dist(t(assay(rld)))
mat <- as.matrix(dists)
rownames(mat) <- colnames(mat) <- as.character(colData$sampleName)
hc <- hclust(dists)

hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(225)
heatmap(mat, Rowv=as.dendrogram(hc), symm=T, trace='none', col=rev(hmcol), margin=c(13,13), main='Heatmap')
dev.copy2pdf(file = "heatmap.pdf")

### PCA 
Pvars <- rowVars(assay(rld))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(1000,length(Pvars)))]
PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
head(PCA)
dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                    condition = cellType)

(qplot(PC1, PC2, data = dataGG, color =  condition,main = "Analyse à composantes principales"))
dev.copy2pdf(file = "PCA.pdf")

### filtration des outliers
outliers <- as.character(subset(colnames(dds), dataGG$PC1 > 0))
outliers
dds_2 = dds [, !(colnames(dds) %in% outliers)] 
rdl_2 = rlogTransformation(dds_2,blind = TRUE)
distsRL = dist(t(assay(rdl_2)))
mat <- as.matrix(distsRL)
rownames(mat) <-  colData(rdl_2)$cellType
hc <- hclust(distsRL)

### heatmap sans outlier
hmcol <- colorRampPalette(brewer.pal(8, "OrRd"))(255)
heatmap(mat, Rowv=as.dendrogram(hc), symm=T, trace='none', col=rev(hmcol), margin=c(13,13), main='Heatmap')
dev.copy2pdf(file = "heatmap_sans_outlier.pdf")

### ACP sans outlier
Pvars <- rowVars(assay(rdl_2))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(1000, length(Pvars)))]

PCA <- prcomp(t(assay(rdl_2)[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

dataGG_w_o = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], PC3 = PCA$x[,3], PC4 = PCA$x[,4], condition = dds_2$cellType)

(qplot(PC1, PC2, data = dataGG_w_o, color =  condition,main = "Analyse à composantes principales sans les outliers"))

dev.copy2pdf(file = "PCA_sans_outlier.pdf")
####  Analyse d'expression différentielle
## graphique de dispersion des gènes 
dds_2 <- estimateDispersions(dds_2)
plotDispEsts(dds_2, main = "Dispersion des gènes sans outlier")
dev.copy2pdf(file = "dispersion_genes_without_outliers.pdf")

####Test statistique de l'expression différentielle
dds_2 = nbinomWaldTest(dds_2)
dds_2res <- results(dds_2, pAdjustMethod = "BH")
table(dds_2res$padj < 0.1)

## Contrôle et correction de p-values
hist(dds_2res$pvalue, col = "lavender", main = "wild vs Mutant", xlab = "p-values")
dds_2res = dds_2res[ !is.na(dds_2res$padj), ]
dds_2res = dds_2res[ !is.na(dds_2res$pvalue), ]
dds_2res <- dds_2res[, -which(names(dds_2res) == "padj")]
fdr_dds_2res = fdrtool(dds_2res$stat, statistic= "normal", plot = T)
fdr_dds_2res$param[1, "sd"]
dds_2res[,"padj"]  <- p.adjust(fdr_dds_2res$pval, method = "BH")

#S'assurer que histogramme des pvalues pour voir leur distribution et si ils ont une distribution uniforme ave
hist(fdr_dds_2res$pval, col = "lavender", main = "wild vs Mutant", xlab = "p-values")
dev.copy2pdf(file = "histogramme_pval.pdf")

#### Visualisation des gènes exprimés
## plotMA
table(dds_2res$padj < 0.1)
plotMA(dds_2res)


##  Volcano plot sur les genes etudiées en biblio 
#par(mfrow = c(1,2))
result_analyse<-dds_2res[c("ENSG00000088256","ENSG00000101019","ENSG00000114770","ENSG00000115524","ENSG00000156052","ENSG00000163930","ENSG00000245694","ENSG00000148848"),]
with(result_analyse, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot on studied genes", xlim=c(-10,10)))
with(subset(dds_2res, padj<0.05 & log2FoldChange>2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(dds_2res, padj<0.05 & log2FoldChange<(-2)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
abline(h = a) 
abline(v = (-2))
abline(v = (2))
dev.copy2pdf(file = "volcano_plot_genes_studied.pdf")

# Volcano plot sur tous les gènes exprimés
with(dds_2res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot for all genes ", xlim=c(-10,10)))
with(subset(dds_2res, padj<0.05 & log2FoldChange>2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(dds_2res, padj<0.05 & log2FoldChange<(-2)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
abline(h = a)
abline(v = (-2))
abline(v = (2))
dev.copy2pdf(file = "volcano_plot_all_genes.pdf")



