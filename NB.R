library("DESeq2")

counts_fp <- "/Users/willtrim/Desktop/RNAseq Bed Rest Data/Bed Rest RNAseq Counts Data.csv"
samples_fp <- "/Users/willtrim/Desktop/RNAseq Bed Rest Data/samples.csv"


setwd("/Users/willtrim/Desktop/RNAseq Bed Rest Data/outputs/deseq2_results")

cts <- as.matrix(read.csv(counts_fp, sep=",",row.names="Geneid"))
samples <- read.csv(samples_fp, sep=",")

# remove an outlier-looking sample (with low RIN number)
samples_ <- samples[!(samples$patient %in% c("I2")),]

samples_$patient <- factor(samples_$patient)
samples_$group <- factor(samples_$group, levels = c("Control", "Cocktail"))
samples_$tmp <- factor(samples_$tmp, levels = c("Pre", "Post"))


cts_ <- cts[, samples_$name]

dds <- DESeqDataSetFromMatrix(countData = cts_,
                              colData = samples_,
                              design = ~ patient + GD.mg.kg.min + tmp)

design_matrix <- model.matrix(~ group + tmp + patient  + tmp:group, data = samples_)
print(design_matrix)
design_matrix <- design_matrix[,-which(colnames(design_matrix) %in% c("patientA2"))]
library(limma)
is.fullrank(design_matrix)
# model.matrix(design(dds), data=as.data.frame(colData(dds)))
dds <- DESeqDataSetFromMatrix(countData = cts_,
                              colData = samples_,
                              design = design_matrix)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))


loggeomeans <- rowMeans(log(counts(dds)))
hist(log(counts(dds)[,1]) - loggeomeans, 
     col="grey", main="", xlab="", breaks=40)

log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
log.norm <- normTransform(dds)

rs <- rowSums(counts(dds))
mypar(1,2)
boxplot(log2(counts(dds)[rs > 0,] + 1)) # not normalized
boxplot(log.norm.counts[rs > 0,]) # normalized
plot(log.norm.counts[,1:2], cex=.1)

rld <- rlog(dds, blind=FALSE)
plot(assay(rld)[,1:2], cex=.1)

library(vsn)
meanSdPlot(log.norm.counts, ranks=FALSE) 
meanSdPlot(assay(rld), ranks=FALSE)

dds <- DESeq(dds)

resultsNames(dds)

res <- results(dds, contrast=c("tmp","Post","Pre"), alpha=0.05)
# res <- results(dds, contrast=list("groupControl.tmpPost", "groupCocktail.tmpPost"), alpha=0.05)
# res <- results(dds, name="groupControl.tmpPost", alpha=0.05)
# res <- results(dds, name="groupCocktail.tmpPost", alpha=0.05)
# res <- results(dds, name="groupCocktail", alpha=0.05)
# res <- results(dds, name="tmpPost", alpha=0.05)
summary(res)


plotPCA(log.norm, intgroup="tmp")
plotPCA(rld, intgroup="tmp")
plotPCA(rld, intgroup="patient")
plotPCA(rld, intgroup="group")

pcaData <- plotPCA(rld, intgroup=c("tmp", "group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=tmp, shape=group.1)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# The MA-plot provides a global view of the differential genes,
# with the log2 fold change on the y-axis over the mean of normalized counts:
plotMA(res, ylim=c(-4,4))

# A p-value histogram:
hist(res$pvalue[res$baseMean > 1], 
     col="grey", border="white", xlab="", ylab="", main="")

# A sorted results table:
resSort <- res[order(res$padj),]
head(resSort)

# Examine the counts for the top gene, sorting by p-value:
gene <- rownames(resSort)[1]
# gene <- topgenes[startsWith(topgenes, "ENSG00000174697")]  # LPT
# gene <- topgenes[startsWith(topgenes, "ENSG00000186185")]  # KIF18B
# gene <- topgenes[startsWith(topgenes, "ENSG00000126787")]  # DLGAP5
# gene <- topgenes[startsWith(topgenes, "ENSG00000174348")]  # PODN
# gene <- topgenes[startsWith(topgenes, "ENSG00000126787")]  # DLGAP5
plotCounts(dds, gene=gene, intgroup=c("tmp", "group"))
plotCounts(dds, gene=gene, intgroup="group")


library(ggplot2)
# gene <- rownames(resSort)[1]
# gene <- "ENSG00000171864.4"
data <- plotCounts(dds, gene=gene, intgroup=c("tmp","group", "patient"), returnData=TRUE, normalized=TRUE)
data$tmp <- factor(data$tmp, levels = c("Pre", "Post"))
data$group <- factor(data$group)

ggplot(data, aes(x=tmp, y=count, col=patient)) +
  geom_point(position=position_jitter(width=.1,height=0)) +
  scale_y_log10()

# Connecting by lines shows the differences which are actually being tested by results given that our0 design includes cell + dex
ggplot(data, aes(x=tmp, y=count, col=group, group=patient)) +
  geom_point() + scale_y_log10() + geom_line()

library(dplyr)
mean_df <- data %>%
  group_by(group) %>%
  summarise(mean_value = mean(count)) # Keep x_var for plotting line

ggplot(data, aes(x=group, y=count, col=group, group=group)) +
  geom_point() + scale_y_log10() + 
  geom_line(data = mean_df, aes(x = group, y = mean_value, color = group), size = 1)
# A heatmap of the top genes:
  
library(pheatmap)
topgenes <- head(rownames(resSort),100)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("tmp","patient")])
pheatmap(mat, annotation_col=df)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

topgenes_ <- sapply(strsplit(topgenes, "\\."), function(x) x[1])

anno <- select(org.Hs.eg.db, keys=topgenes_,
               columns=c("SYMBOL","GENENAME"), 
               keytype="ENSEMBL")

# to csv
write.csv(res, "deseq2_results.csv", row.names = TRUE)
# write.csv(res, "deseq2_results_groupCocktail_tmpPost.csv", row.names = TRUE)

normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, "normalized_counts.csv", row.names = TRUE)


# access to values

mu <- assays(dds)[["mu"]]
head(mu)
mu[gene,]

coefs <- coef(dds)
head(coefs)
coefs[gene,]
