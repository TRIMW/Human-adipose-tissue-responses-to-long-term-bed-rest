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


deseq2_results_fp <- "/n/groups/kirschner/Will/BedRest/deseq2/star_rsem_tmp_post_pre.csv"
deseq2_results <- read.csv(deseq2_results_fp)
deseq2_results$ensembl_gene_id <- sapply(strsplit(deseq2_results$X, "\\."), head, 1)
deseq2_results <- merge(genes_unique, deseq2_results, by = "ensembl_gene_id")

deseq2_ranks <- deseq2_results$stat
names(deseq2_ranks) <- deseq2_results$entrezgene_id
deseq2_ranks <- deseq2_ranks[!is.na(deseq2_ranks)]


# msigdb.hs = getMsigdb(org = 'hs', id = 'EZID', version = '7.5.1')
# msigdb.hs = appendKEGG(msigdb.hs)

msigdb_gmt_fp <- "/n/groups/kirschner/reference/msigdb/msigdb_v2025.1.Hs_files_to_download_locally/msigdb_v2025.1.Hs_GMTs/"


### Reactome

set.seed(42)

gmt.file <- paste(msigdb_gmt_fp, "c2.cp.reactome.v2025.1.Hs.entrez.gmt", sep="")
pathways <- gmtPathways(gmt.file)
str(head(pathways))

fgseaRes_reactome <- fgsea(pathways = pathways, 
                  stats    = deseq2_ranks,
                  minSize  = 10,
                  maxSize  = 500)

head(fgseaRes_reactome[order(pval), ])
sum(fgseaRes_reactome[, padj < 0.01])

# Plot the most significantly enriched pathway
plotEnrichment(pathways[[head(fgseaRes_reactome[order(pval), ], 1)$pathway]],
               deseq2_ranks) + labs(title=head(fgseaRes_reactome[order(pval), ], 1)$pathway)


pathway <- "REACTOME_PENTOSE_PHOSPHATE_PATHWAY"
plotEnrichment(pathways[[pathway]],
               deseq2_ranks) + labs(title=pathway)

# Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively.
topPathwaysUp <- fgseaRes_reactome[ES > 0][head(order(pval), n=25), pathway]
topPathwaysDown <- fgseaRes_reactome[ES < 0][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], deseq2_ranks, fgseaRes_reactome, 
              gseaParam = 0.5)

# fgsea_reactome[order(pval),][1,]$leadingEdge
## 

###  Hallmark 

set.seed(42)

hallmark.gmt.file <- paste(msigdb_gmt_fp, "h.all.v2025.1.Hs.entrez.gmt", sep="")
hallmark.pathways <- gmtPathways(hallmark.gmt.file)
str(head(hallmark.pathways))

fgseaRes_hallmark <- fgsea(pathways = hallmark.pathways, 
                           stats    = deseq2_ranks,
                           minSize  = 15,
                           maxSize  = 500)

head(fgseaRes_hallmark[order(pval), ])
sum(fgseaRes_hallmark[, padj < 0.01])

# Plot the most significantly enriched pathway
plotEnrichment(hallmark.pathways[[head(fgseaRes_hallmark[order(pval), ], 1)$pathway]],
               deseq2_ranks) + labs(title=head(fgseaRes_hallmark[order(pval), ], 1)$pathway)

# Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively.
topPathwaysUp <- fgseaRes_hallmark[ES > 0][head(order(pval), n=15), pathway]
topPathwaysDown <- fgseaRes_hallmark[ES < 0][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(hallmark.pathways[topPathways], deseq2_ranks, fgseaRes_hallmark, 
              gseaParam = 0.5)

# fgseaRes_hallmark[order(pval),][1,]$leadingEdge

### Wikipathways
wiki.gmt.file <- paste(msigdb_gmt_fp, "c2.cp.wikipathways.v2025.1.Hs.entrez.gmt", sep="")
wiki.pathways <- gmtPathways(wiki.gmt.file)
str(head(wiki.pathways))

fgseaRes_wiki <- fgsea(pathways = wiki.pathways, 
                           stats    = deseq2_ranks,
                           minSize  = 15,
                           maxSize  = 500)

head(fgseaRes_wiki[order(pval), ])
sum(fgseaRes_wiki[, padj < 0.01])

# Plot the most significantly enriched pathway
plotEnrichment(wiki.pathways[[head(fgseaRes_wiki[order(pval), ], 1)$pathway]],
               deseq2_ranks) + labs(title=head(fgseaRes_wiki[order(pval), ], 1)$pathway)

# Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively.
topPathwaysUp <- fgseaRes_wiki[ES > 0][head(order(pval), n=15), pathway]
topPathwaysDown <- fgseaRes_wiki[ES < 0][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(wiki.pathways[topPathways], deseq2_ranks, fgseaRes_wiki, 
              gseaParam = 0.5)

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

fgseaRes_kegg <- fgsea(
  pathways = kegg_pathways, 
  stats    = deseq2_ranks,
  minSize  = 15,
  maxSize  = 500,
  eps=0
)

head(fgseaRes_kegg[order(pval), ])
sum(fgseaRes_kegg[, padj < 0.01])

# Plot the most significantly enriched pathway
plotEnrichment(kegg_pathways[[head(fgseaRes_kegg[order(pval), ], 1)$pathway]],
               deseq2_ranks) + labs(title=head(fgseaRes_kegg[order(pval), ], 1)$pathway)

# Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively.
topPathwaysUp <- fgseaRes_kegg[ES > 0][head(order(pval), n=15), pathway]
topPathwaysDown <- fgseaRes_kegg[ES < 0][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(kegg_pathways[topPathways], deseq2_ranks, fgseaRes_kegg, 
              gseaParam = 0.5)

# fgseaRes_hallmark[order(pval),][1,]$leadingEdge

### GTRD

gtrd.gmt.file <- paste(msigdb_gmt_fp, "c3.tft.gtrd.v2025.1.Hs.entrez.gmt", sep="")
gtrd.pathways <- gmtPathways(gtrd.gmt.file)
str(head(gtrd.pathways))

fgseaRes_gtrd <- fgsea(pathways = gtrd.pathways, 
                       stats    = deseq2_ranks,
                       minSize  = 15,
                       maxSize  = 500)

head(fgseaRes_gtrd[order(pval), ])
sum(fgseaRes_gtrd[, padj < 0.01])

# Plot the most significantly enriched pathway
plotEnrichment(gtrd.pathways[[head(fgseaRes_gtrd[order(pval), ], 1)$pathway]],
               deseq2_ranks) + labs(title=head(fgseaRes_gtrd[order(pval), ], 1)$pathway)

# Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively.
topPathwaysUp <- fgseaRes_gtrd[ES > 0][head(order(pval), n=5), pathway]
topPathwaysDown <- fgseaRes_gtrd[ES < 0][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(gtrd.pathways[topPathways], deseq2_ranks, fgseaRes_gtrd, 
              gseaParam = 0.5)


### GO

#### BP
bp.gmt.file <- paste(msigdb_gmt_fp, "c5.go.bp.v2025.1.Hs.entrez.gmt", sep="")
bp.pathways <- gmtPathways(bp.gmt.file)
str(head(bp.pathways))

fgseaRes_bp <- fgsea(pathways = bp.pathways, 
                       stats    = deseq2_ranks,
                       minSize  = 15,
                       maxSize  = 500)

head(fgseaRes_bp[order(pval), ])
sum(fgseaRes_bp[, padj < 0.01])

# Plot the most significantly enriched pathway
plotEnrichment(bp.pathways[[head(fgseaRes_bp[order(pval), ], 1)$pathway]],
               deseq2_ranks) + labs(title=head(fgseaRes_bp[order(pval), ], 1)$pathway)

# Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively.
topPathwaysUp <- fgseaRes_bp[ES > 0][head(order(pval), n=15), pathway]
topPathwaysDown <- fgseaRes_bp[ES < 0][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(bp.pathways[topPathways], deseq2_ranks, fgseaRes_bp, 
              gseaParam = 0.5)

#### CC
cc.gmt.file <- paste(msigdb_gmt_fp, "c5.go.cc.v2025.1.Hs.entrez.gmt", sep="")
cc.pathways <- gmtPathways(cc.gmt.file)
str(head(cc.pathways))

fgseaRes_cc <- fgsea(pathways = cc.pathways, 
                       stats    = deseq2_ranks,
                       minSize  = 15,
                       maxSize  = 500)

head(fgseaRes_cc[order(pval), ])
sum(fgseaRes_cc[, padj < 0.01])

# Plot the most significantly enriched pathway
plotEnrichment(cc.pathways[[head(fgseaRes_cc[order(pval), ], 1)$pathway]],
               deseq2_ranks) + labs(title=head(fgseaRes_cc[order(pval), ], 1)$pathway)

# Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively.
topPathwaysUp <- fgseaRes_cc[ES > 0][head(order(pval), n=15), pathway]
topPathwaysDown <- fgseaRes_cc[ES < 0][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(cc.pathways[topPathways], deseq2_ranks, fgseaRes_cc, 
              gseaParam = 0.5)

#### MF

mf.gmt.file <- paste(msigdb_gmt_fp, "c5.go.mf.v2025.1.Hs.entrez.gmt", sep="")
mf.pathways <- gmtPathways(mf.gmt.file)
str(head(mf.pathways))

fgseaRes_mf <- fgsea(pathways = mf.pathways, 
                       stats    = deseq2_ranks,
                       minSize  = 15,
                       maxSize  = 500)

head(fgseaRes_mf[order(pval), ])
sum(fgseaRes_mf[, padj < 0.01])

# Plot the most significantly enriched pathway
plotEnrichment(mf.pathways[[head(fgseaRes_mf[order(pval), ], 1)$pathway]],
               deseq2_ranks) + labs(title=head(fgseaRes_mf[order(pval), ], 1)$pathway)

# Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively.
topPathwaysUp <- fgseaRes_mf[ES > 0][head(order(pval), n=15), pathway]
topPathwaysDown <- fgseaRes_mf[ES < 0][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(mf.pathways[topPathways], deseq2_ranks, fgseaRes_mf, 
              gseaParam = 0.5)
