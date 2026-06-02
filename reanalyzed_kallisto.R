# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("tximport")
# BiocManager::install("rhdf5")
# BiocManager::install("limma")

library(tximport)
library(DESeq2)

samples_fp <- '/n/groups/kirschner/Will/BedRest/deseq2/samples_kallisto_tximport.csv'
samples <- read.csv(samples_fp, sep=",")
files <- samples$filepath
names(files) <- samples$sample_name
rownames(samples) <- samples$sample_name

# remove an outlier-looking sample (with low RIN number)
# samples <- samples[!(samples$patient %in% c("I2")),]

samples$patient <- factor(samples$patient)
samples$group <- factor(samples$group, levels = c("Control", "Cocktail"))
samples$tmp <- factor(samples$tmp, levels = c("Pre", "Post"))

tx2gene_fp <- '/n/groups/kirschner/Will/BedRest/nfcore_rnaseq_3.14.0/kallisto/tx2gene.tsv'
tx2gene <- read.csv(tx2gene_fp, sep="\t", header=FALSE)

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)


# design_matrix <- model.matrix(~ group + patient + tmp:group, data = samples)
# print(design_matrix)
# design_matrix <- design_matrix[,-which(colnames(design_matrix) %in% c("patientA2"))]
# library(limma)
# is.fullrank(design_matrix)
# model.matrix(design(dds), data=as.data.frame(colData(dds)))

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ patient + tmp)
# dds <- DESeqDataSetFromTximport(txi, colData = samples, design = design_matrix)
dds <- DESeq(dds)

resultsNames(dds)

res <- results(dds, contrast=c("tmp","Post","Pre"), alpha=0.05)
# res <- results(dds, contrast=list("groupCocktail.tmpPost", "groupControl.tmpPost"), alpha=0.05)
# res <- results(dds, name="groupControl.tmpPost", alpha=0.05)
summary(res)

# A p-value histogram:
hist(res$pvalue[res$baseMean > 1], 
     col="grey", border="white", xlab="", ylab="", main="")

# A sorted results table:
resSort <- res[order(res$padj),]
head(resSort)
gene <- rownames(resSort)[5]
plotCounts(dds, gene=gene, intgroup="tmp")

library(ggplot2)
data <- plotCounts(dds, gene=gene, intgroup=c("tmp", "patient", "group"), returnData=TRUE, normalized=TRUE)
# data$tmp <- factor(data$tmp, levels = c("Pre", "Post"))
# data$group <- factor(data$group)

ggplot(data, aes(x=tmp, y=count, col=patient)) +
  geom_point(position=position_jitter(width=.1,height=0)) +
  scale_y_log10()

# Connecting by lines shows the differences which are actually being tested by results given that our0 design includes cell + dex
ggplot(data, aes(x=tmp, y=count, col=patient, group=patient)) +
  geom_point() + scale_y_log10() + geom_line()




setwd("/n/groups/kirschner/Will/BedRest/deseq2")

write.csv(res, "deseq2_results_kallisto_tmp_post_pre.csv", row.names = TRUE)

normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, "normalized_counts_kallisto_tmp_post_pre.csv", row.names = TRUE)

