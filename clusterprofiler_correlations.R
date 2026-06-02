library(biomaRt)
library(fgsea)
library(msigdb)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
devtools::load_all("~/projs/enrichplot")

counts_fp <- "/n/groups/kirschner/Will/BedRest/nfcore_rnaseq_3.14.0/star_rsem/rsem.merged.gene_counts.tsv"
counts <- read.csv(counts_fp, sep="\t")
counts$transcript_id.s. <- NULL


gene_ids <- sapply(strsplit(counts$gene_id, "\\."), head, 1)

# Get the corresponding Entrez IDs
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id", "entrezgene_id"),
               values = gene_ids,
               mart = mart)

genes_unique <- genes %>% distinct(ensembl_gene_id, .keep_all = TRUE)

genes_unique <- genes_unique %>%
  group_by(entrezgene_id) %>%
  filter(n() == 1) %>%
  ungroup()

counts$ensembl_gene_id <- gene_ids
nrow(counts)
counts <- merge(genes_unique, counts, by = "ensembl_gene_id")
nrow(counts)
rownames(counts) <- counts$entrezgene_id
counts$ensembl_gene_id <- NULL
counts$entrezgene_id <- NULL
counts$gene_id <- NULL

##### Correlations input

results_fps <- list.files(
  path = "/n/groups/kirschner/Will/BedRest/deseq2/cor", 
  pattern = ".*.csv$",
  full.names = TRUE # returns the full path, not just the filename
)

pheno_ranks <- list()

for (results_fp in results_fps) { 
  results <- read.csv(results_fp)
  results$ensembl_gene_id <- sapply(strsplit(results$gene, "\\."), head, 1)
  results <- merge(genes_unique, results, by = "ensembl_gene_id")
  
  pheno <- strsplit(basename(results_fp), "\\.")[[1]][[1]]
  
  results <- results[order(-results$r),]
  ranks <- results$r
  names(ranks) <- results$entrezgene_id
  ranks <- ranks[!is.na(ranks)]
  
  pheno_ranks[[pheno]] <- ranks
}


all_phenotypes <- c(
  "TAG",
  "FFA",
  "glucose_art",
  "glucose_ven",
  "insulin_art",
  "insulin_ven",
  "cho_ox",
  "fat_ox"
)

pheno_tmp <- unlist(sapply(all_phenotypes, function(x) {return(c(paste0("pre_", x), paste0("post_", x)))}, simplify=FALSE))

pheno_ranks <- pheno_ranks[pheno_tmp]

## KEGG
xx <- compareCluster(pheno_ranks,
                     fun="gseKEGG",
                     organism="hsa",
                     minGSSize    = 15,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.05,
                     verbose      = TRUE,
                     eps=0,
                     nPermSimple = 100000
)


xx@compareClusterResult$ID_ <- xx@compareClusterResult$ID
xx@compareClusterResult$ID <- xx@compareClusterResult$Description

xx <- pairwise_termsim(xx)       
library(stringr)
xx@compareClusterResult$setSize <- as.numeric(str_match(xx@compareClusterResult$leading_edge, "tags=([0-9]+)%")[,2]) / 100

# remove_sets <- c("mitochondrion organization", "")
# xx@compareClusterResult <- filter(xx@compareClusterResult, Description %in% remove_sets)

treeplot(xx,
         showCategory = 10,
         #size_var = c("Count"), # "setSize"),
         cluster_panel="dotplot",
         cluster_method = "complete",
         color="NES",
         nCluster=7,
         fontsize_tiplab = 4,
         fontsize_cladelab = 4,
         group_color = NULL,
         hexpand=0.05,
         tiplab_offset=40,
         cladelab_offset=70,
         colnames_angle=45
         # nWords=0,
         # cluster.params = list(method = "ward.D", n = 5, color = "NES", label_words_n = 0,
         # label_format = 30),
         # hilight.params = list(hilight = TRUE, align = "both"),
         # clusterPanel.params = list(clusterPanel = "heatMap", pie = "equal", legend_n = 3,
         # colnames_angle = 0),
         # offset.params = list(bar_tree = 1, cladelab=10, tiplab = 25, extend = 1, hexpand = 2),
) + 
  scale_size_continuous(range = c(1, 4), name="LE %") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0)))

