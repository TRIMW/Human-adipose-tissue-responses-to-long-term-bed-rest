# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("tximport")
# BiocManager::install("rhdf5")
# BiocManager::install("limma")

library(tximport)
library(DESeq2)


samples_fp <- '/n/groups/kirschner/Will/BedRest/deseq2/samples_star_rsem_tximport.csv'
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

txi <- tximport(files, type = "rsem", tx2gene = tx2gene)


design_matrix <- model.matrix(~ group + patient  + tmp:group, data = samples)
print(design_matrix)
# design_matrix <- design_matrix[,-which(colnames(design_matrix) %in% c("patientA2"))]
design_matrix <- design_matrix[,-which(colnames(design_matrix) %in% c("groupCocktail"))]
library(limma)
is.fullrank(design_matrix)
# model.matrix(design(dds), data=as.data.frame(colData(dds)))

# dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ patient + tmp)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = design_matrix)

# filter low counts
keep <- rowSums(counts(dds) >= 5) >= 3
dds <- dds[keep,]
dds <- DESeq(dds)

resultsNames(dds)

# res <- results(dds, contrast=c("tmp","Post","Pre"), alpha=0.05)
res <- results(dds, contrast=list("groupCocktail.tmpPost", "groupControl.tmpPost"), alpha=0.05)
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

write.csv(res, "star_rsem_tmp_post_pre.csv", row.names = TRUE)

normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, "normalized_counts_star_rsem_tmp_post_pre.csv", row.names = TRUE)

rld <- rlog(dds, blind = TRUE) # Set blind=TRUE for unbiased QA
rlog_matrix <- assay(rld)
write.csv(rlog_matrix, "rlog_counts_star_rsem.csv", row.names = TRUE)

# PCA
pcaData <- plotPCA(rld, intgroup=c("tmp"), returnData=TRUE, ntop=1500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$participant <- sapply(strsplit(pcaData$name, "_"), head, 1)
pcaData_pre <- pcaData[pcaData$tmp == "Pre",]
pcaData_post <- pcaData[pcaData$tmp == "Post",]
pcaData_pre_post <- merge(pcaData_pre, pcaData_post, by = "participant", suffixes = c("_pre", "_post"))

library(ggplot2)

p <- ggplot(pcaData, aes(PC1, PC2, color=tmp)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio = 2.2) +
  stat_ellipse(geom = "polygon",
               aes(fill = tmp, group=tmp), 
               alpha = 0.1) +
  scale_color_manual(values = c("Pre" = "#999797", "Post" = "#a73bf5")) +
  scale_fill_manual(values = c("Pre" = "#ACACAC", "Post" = "#BE6FF7")) +
  theme_minimal(base_size=38) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 36), #, face = "bold"), # Change legend title font size and style
        legend.text = element_text(size = 32) # Change legend text (labels) font size
  ) + guides(fill = guide_legend(override.aes = list(size = 6, shape = 15, linetype = 0, fill=NA))) +
  geom_segment(
    data=pcaData_pre_post, aes(x = PC1_pre, y = PC2_pre, xend = PC1_post, yend = PC2_post),
    arrow = arrow(length = unit(0.2, "cm")), # Customize the arrow head length
    color = "black",                          # Optional: Set the color
    size = 0.5                                 # Optional: Set the line thickness
  )
p

ggplot(pcaData) +
  geom_violin(aes(y=PC1, x=tmp))

ggplot(pcaData) +
  geom_violin(aes(y=PC2, x=tmp))


# Volcano plot
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)

res_ids <- sapply(strsplit(rownames(resSort), "\\."), head, 1)

res_names <- bitr(res_ids, fromType="ENSEMBL", 
                  toType="SYMBOL", 
                  OrgDb="org.Hs.eg.db") #NA values drop by default

res_names <- res_names %>%
  distinct(ENSEMBL, .keep_all = TRUE)

resSort[(!is.na(resSort$pvalue))&(resSort$log2FoldChange<0)&(resSort$pvalue<0.01),]
resSort[(!is.na(resSort$pvalue))&(resSort$log2FoldChange>0)&(resSort$pvalue<0.01),]

rownames(res_names) <- res_names$ENSEMBL

resSort$symbol <- res_names[res_ids, "SYMBOL"]

pval_threshold <- 0.05
keyvals <- ifelse(
  (resSort$padj < pval_threshold) & (resSort$log2FoldChange > 0), 'red',
  ifelse((resSort$padj < pval_threshold) & (resSort$log2FoldChange < 0), 'royalblue',
         '#4d4c4c'))
keyvals[is.na(keyvals)] <- '#4d4c4c'
names(keyvals)[keyvals == 'red'] <- paste('Upregulated (Padj < ', as.character(pval_threshold), ')')
names(keyvals)[keyvals == '#4d4c4c'] <- 'Non-DEGs'
names(keyvals)[keyvals == 'royalblue'] <- paste('Downregulated (Padj < ', as.character(pval_threshold), ')')

resSort_ <- resSort[(!is.na(resSort$padj)) & (resSort$padj < pval_threshold),] #(abs(resSort$log2FoldChange) > 0.5) & 


p <- EnhancedVolcano(
  resSort,
  lab=resSort$symbol,
  x="log2FoldChange",
  y="padj",
  pCutoff=0.05,
  FCcutoff=0,
  title="",
  subtitle="",
  drawConnectors=TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  arrowheads = FALSE,
  pointSize = 1.2,
  labSize=5,
  selectLab=resSort_$symbol[1:50],
  # colAlpha = 1,
  colCustom = keyvals,
  cutoffLineType = 'blank',
  ylim = c(0, 22),
  xlim = c(-4, 4),
  hline = c(0.05),
  hlineCol = c('#333232'),
  hlineType = c('longdash'),
  hlineWidth = c(0.5),
  # gridlines.major = FALSE,
  gridlines.minor = FALSE,
  border = 'partial',
  borderWidth = 1,
  borderColour = 'black',
  legendPosition = 'right',
  legendLabSize = 15,
  legendIconSize = 5.0,
  caption=NULL,
  # min.segment.length=0.5
) + 
  labs(color = NULL) +
  theme_minimal(base_size=28) +
  theme(legend.title = element_text(size = 26), #, face = "bold"), # Change legend title font size and style
        legend.text = element_text(size = 18) # Change legend text (labels) font size
  ) #+ guides(fill = guide_legend(override.aes = list(size = 6, shape = 15, linetype = 0, fill=NA)))


p
setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2")
ggsave("volcano.png", plot = p, dpi = 400)