# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("tximport")

library(tximport)

# Example: Adjust paths to your Kallisto output directories
samples <- data.frame(
  sample_name = c("sample1", "sample2", "sample3"),
  path = c("path/to/sample1/abundance.h5", "path/to/sample2/abundance.h5", "path/to/sample3/abundance.h5")
)
files <- samples$path
names(files) <- samples$sample_name

# Example using GenomicFeatures (adjust for your specific organism and GTF)
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("path/to/your/annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)

# Example with DESeq2
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
dds <- DESeq(dds)
