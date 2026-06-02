# library(devtools)
# install_github("ctlab/fgsea")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")
# BiocManager::install("msigdb")
# BiocManager::install("EnrichmentBrowser")
# BiocManager::install("clusterProfiler")


library(biomaRt)
library(fgsea)
library(msigdb)
library(dplyr)
library(ggplot2)

counts_fp <- "/Users/willtrim/Documents/projs/bedrest/outputs/rsem.merged.gene_counts.tsv"
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
  path = "/Users/willtrim/Documents/projs/bedrest/outputs/delta_cors", 
  pattern = ".*.csv$",
  full.names = TRUE # returns the full path, not just the filename
)

pheno_ranks <- list()

for (results_fp in results_fps) { 
  results <- read.csv(results_fp)
  results$ensembl_gene_id <- sapply(strsplit(results$gene, "\\."), head, 1)
  results <- merge(genes_unique, results, by = "ensembl_gene_id")
  
  pheno <- strsplit(basename(results_fp), "\\.")[[1]][[1]]
  
  ranks <- results$r
  names(ranks) <- results$entrezgene_id
  ranks <- ranks[!is.na(ranks)]
  
  pheno_ranks[[pheno]] <- ranks
}

# msigdb.hs = getMsigdb(org = 'hs', id = 'EZID', version = '7.5.1')
# msigdb.hs = appendKEGG(msigdb.hs)

msigdb_gmt_fp <- "/Users/willtrim/Documents/reference/msigdb/msigdb_v2025.1.Hs_GMTs/"
min_size = 10
max_size = 500
nPermSimple = 100000

### Reactome

set.seed(42)

gmt.file <- paste(msigdb_gmt_fp, "c2.cp.reactome.v2025.1.Hs.entrez.gmt", sep="")
react.pathways <- gmtPathways(gmt.file)
str(head(react.pathways))

react_res <- list()
pathways <- react.pathways

for (pheno_name in names(pheno_ranks)) {
  fgseaRes_reactome <- fgsea(pathways = pathways, 
                             stats    = pheno_ranks[[pheno_name]],
                             minSize  = min_size,
                             maxSize  = max_size,
                             eps=0,
                             nPermSimple = nPermSimple
                             )
  print(pheno_name)
  print(head(fgseaRes_reactome[order(pval), ]))
  # print(sum(fgseaRes_reactome[, padj < 0.01]))
  
  react_res[[pheno_name]] <- fgseaRes_reactome
}


## 

###  Hallmark 

set.seed(42)

hallmark.gmt.file <- paste(msigdb_gmt_fp, "h.all.v2025.1.Hs.entrez.gmt", sep="")
hallmark.pathways <- gmtPathways(hallmark.gmt.file)
str(head(hallmark.pathways))

hallmark_res <- list()
pathways <- hallmark.pathways

for (pheno_name in names(pheno_ranks)) {
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = pheno_ranks[[pheno_name]],
                    minSize  = min_size,
                    maxSize  = max_size,
                    eps=0,
                    nPermSimple = nPermSimple
                    )
  print(pheno_name)
  print(head(fgseaRes[order(pval), ]))
  
  hallmark_res[[pheno_name]] <- fgseaRes
}


### Wikipathways
wiki.gmt.file <- paste(msigdb_gmt_fp, "c2.cp.wikipathways.v2025.1.Hs.entrez.gmt", sep="")
wiki.pathways <- gmtPathways(wiki.gmt.file)
str(head(wiki.pathways))

wiki_res <- list()
pathways <- wiki.pathways

for (pheno_name in names(pheno_ranks)) {
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = pheno_ranks[[pheno_name]],
                    minSize  = min_size,
                    maxSize  = max_size,
                    eps=0,
                    nPermSimple = nPermSimple
                    )
  print(pheno_name)
  print(head(fgseaRes[order(pval), ]))
  wiki_res[[pheno_name]] <- fgseaRes
}

### KEGG

# library(KEGGREST)
# listDatabases()
# queryables <- c(listDatabases(), org[,1], org[,2])
# kegg_pathways <- keggLink("pathway", "hsa")

# library(EnrichmentBrowser)
library(clusterProfiler)

kegg_pathways_df <- download_KEGG(species = "hsa", keggType = "KEGG", keyType = "kegg")

kegg_pathways_list <- kegg_pathways_df["KEGGPATHID2EXTID"][[1]] %>%
  group_by(from) %>%
  group_split()

library(purrr)
kegg_pathways <- map(kegg_pathways_list, ~pull(.x, to))
kegg_pathway_ids <-  map(kegg_pathways_list, ~ .x[[1, 1]])
kegg_pathways_names_map <- kegg_pathways_df[["KEGGPATHID2NAME"]]$to
names(kegg_pathways_names_map) <- kegg_pathways_df[["KEGGPATHID2NAME"]]$from
kegg_pathways_names <- lapply(kegg_pathway_ids, function(x) {kegg_pathways_names_map[[x]]})
names(kegg_pathways) <- kegg_pathways_names

kegg_res <- list()

for (pheno_name in names(pheno_ranks)) {
  fgseaRes <- fgsea(pathways = kegg_pathways, 
                    stats    = pheno_ranks[[pheno_name]],
                    minSize  = min_size,
                    maxSize  = max_size,
                    eps=0,
                    nPermSimple = nPermSimple
  )
  print(pheno_name)
  print(head(fgseaRes[order(pval), ]))
  kegg_res[[pheno_name]] <- fgseaRes
}

### GTRD

gtrd.gmt.file <- paste(msigdb_gmt_fp, "c3.tft.gtrd.v2025.1.Hs.entrez.gmt", sep="")
gtrd.pathways <- gmtPathways(gtrd.gmt.file)
str(head(gtrd.pathways))

gtrd_res <- list()
pathways <- gtrd.pathways

for (pheno_name in names(pheno_ranks)) {
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = pheno_ranks[[pheno_name]],
                    minSize  = min_size,
                    maxSize  = max_size,
                    eps=0,
                    nPermSimple = nPermSimple)
  print(pheno_name)
  print(head(fgseaRes[order(pval), ]))
  gtrd_res[[pheno_name]] <- fgseaRes
}

### GO

#### BP
bp.gmt.file <- paste(msigdb_gmt_fp, "c5.go.bp.v2025.1.Hs.entrez.gmt", sep="")
bp.pathways <- gmtPathways(bp.gmt.file)
str(head(bp.pathways))

bp_res <- list()
pathways <- bp.pathways

for (pheno_name in names(pheno_ranks)) {
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = pheno_ranks[[pheno_name]],
                    minSize  = min_size,
                    maxSize  = max_size,
                    eps=0,
                    nPermSimple = nPermSimple)
  print(pheno_name)
  print(head(fgseaRes[order(pval), ]))
  bp_res[[pheno_name]] <- fgseaRes
}


#### CC
cc.gmt.file <- paste(msigdb_gmt_fp, "c5.go.cc.v2025.1.Hs.entrez.gmt", sep="")
cc.pathways <- gmtPathways(cc.gmt.file)
str(head(pathways))

cc_res <- list()
pathways <- cc.pathways

for (pheno_name in names(pheno_ranks)) {
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = pheno_ranks[[pheno_name]],
                    minSize  = min_size,
                    maxSize  = max_size,
                    eps=0,
                    nPermSimple = nPermSimple)
  print(pheno_name)
  print(head(fgseaRes[order(pval), ]))
  cc_res[[pheno_name]] <- fgseaRes
}

#### MF

mf.gmt.file <- paste(msigdb_gmt_fp, "c5.go.mf.v2025.1.Hs.entrez.gmt", sep="")
mf.pathways <- gmtPathways(mf.gmt.file)
str(head(mf.pathways))

mf_res <- list()
pathways <- mf.pathways

for (pheno_name in names(pheno_ranks)) {
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = pheno_ranks[[pheno_name]],
                    minSize  = min_size,
                    maxSize  = max_size,
                    eps=0,
                    nPermSimple = nPermSimple)
  print(pheno_name)
  print(head(fgseaRes[order(pval), ]))
  mf_res[[pheno_name]] <- fgseaRes
}


library(ggplot2)
library(cowplot)
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")
library(ggpubr)

setwd("/n/groups/kirschner/Will/BedRest/deseq2/cor/fgsea_figs_small_gs")

all_results <- list(
  # "reactome"=react_res,
  "KEGG"=kegg_res
  # "Hallmark"=hallmark_res,
  # "GOBP"=bp_res,
  # "Wiki"=wiki_res
  # "GTRD"=gtrd_res,
  # "GOCC"=cc_res,
  # "GOMF"=mf_res
  )

all_pathways <- list(
  # "reactome"=react.pathways,
  "KEGG"=kegg_pathways
  # "Hallmark"=hallmark.pathways,
  # "GOBP"=bp.pathways,
  # "Wiki"=wiki.pathways
  # "GTRD"=gtrd.pathways,
  # "GOCC"=cc.pathways,
  # "GOMF"=mf.pathways
)

all_phenotypes <- c(
  "TAG",
  "FFA",
  "glucose_art",
  "glucose_ven",
  "insulin_art",
  "insulin_ven",
  "cho_ox",
  "fat_ox",
  "glucose_diposal"
  )



plot_table <- function(db_name, pheno_name, pval_threshold=0.05) {
  
  pathways <- all_pathways[[db_name]]
  
  
  pre_pheno <- paste("pre_", pheno_name, sep="")
  post_pheno <- paste("post_", pheno_name, sep="")
  
  
  fgseaRes <- all_results[[db_name]][[pre_pheno]]
  topPathwaysUp <- fgseaRes[(ES > 0) & (padj < pval_threshold)][head(order(pval), n=20), pathway]
  topPathwaysDown <- fgseaRes[(ES < 0) & (padj < pval_threshold)][head(order(pval), n=20), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  top <- length(topPathways)
  
  pre <- plotGseaTable(pathways[topPathways], pheno_ranks[[pre_pheno]], fgseaRes, 
                       gseaParam = 0.5)
  
  fgseaRes <- all_results[[db_name]][[post_pheno]]
  topPathwaysUp <- fgseaRes[(ES > 0) & (padj < pval_threshold)][head(order(pval), n=20), pathway]
  topPathwaysDown <- fgseaRes[(ES < 0) & (padj < pval_threshold)][head(order(pval), n=20), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  bottom <- length(topPathways)
  
  post <- plotGseaTable(pathways[topPathways], pheno_ranks[[post_pheno]], fgseaRes, 
                        gseaParam = 0.5)
  
  p <- ggarrange(pre, post,
                 ncol = 1, nrow = 2, heights=c(top, bottom))
  
  print(
    annotate_figure(
      p,
      fig.lab = paste(pheno_name, db_name, sep=": "),
      # fig.lab.face = "bold"
    )
  )
  
}


for (pheno_name in all_phenotypes) {
    
    for (db_name in names(all_pathways)) {
      jpeg(filename = paste(pheno_name, "_", db_name, ".jpg", sep=""),
           width = 1100,
           height = 1200,
           units = "px",
           pointsize = 12,
           quality = 75 # JPEG specific quality setting
      )
      plot_table(db_name, pheno_name)
      dev.off()
      # break
    }
  # break
}


setwd("/Users/willtrim/Documents/projs/bedrest/outputs/delta_cors/fgsea_res")
library(data.table)

for (pheno_name in all_phenotypes) {
  
  for (db_name in names(all_pathways)) {
    fgseaRes <- all_results[[db_name]][[pheno_name]]
    fwrite(
      fgseaRes,
      paste(pheno_name, ".", db_name, ".csv", sep=""),
      row.names = FALSE
    )
  }
}
