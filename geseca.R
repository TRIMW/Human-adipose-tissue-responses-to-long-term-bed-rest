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

counts_fp <- "/n/groups/kirschner/Will/BedRest/deseq2/rlog_counts_star_rsem.csv"
counts <- read.csv(counts_fp, sep=",")
counts <- counts[,!grepl("I2", colnames(counts))]

gene_ids <- sapply(strsplit(counts$X, "\\."), head, 1)

# Get the corresponding Entrez IDs
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id", "entrezgene_id"),
               values = gene_ids,
               mart = mart)

genes <- genes[order(genes$entrezgene_id),]

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
counts$X <- NULL
nrow(counts)
counts <- counts[rowSums(counts) > 10,]
nrow(counts)


colnames(counts) <- lapply(strsplit(colnames(counts), split = "_"), function(x) paste(x[[2]], x[[1]], sep="_"))


msigdb_gmt_fp <- "/n/groups/kirschner/reference/msigdb/msigdb_v2025.1.Hs_files_to_download_locally/msigdb_v2025.1.Hs_GMTs/"


### Reactome


gmt.file <- paste(msigdb_gmt_fp, "c2.cp.reactome.v2025.1.Hs.entrez.gmt", sep="")
pathways <- gmtPathways(gmt.file)
str(head(pathways))

###  Hallmark 

hallmark.gmt.file <- paste(msigdb_gmt_fp, "h.all.v2025.1.Hs.entrez.gmt", sep="")
hallmark.pathways <- gmtPathways(hallmark.gmt.file)
str(head(hallmark.pathways))

### Wikipathways
wiki.gmt.file <- paste(msigdb_gmt_fp, "c2.cp.wikipathways.v2025.1.Hs.entrez.gmt", sep="")
wiki.pathways <- gmtPathways(wiki.gmt.file)
str(head(wiki.pathways))


### KEGG

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

### GTRD

gtrd.gmt.file <- paste(msigdb_gmt_fp, "c3.tft.gtrd.v2025.1.Hs.entrez.gmt", sep="")
gtrd.pathways <- gmtPathways(gtrd.gmt.file)
str(head(gtrd.pathways))


### GO

#### BP
bp.gmt.file <- paste(msigdb_gmt_fp, "c5.go.bp.v2025.1.Hs.entrez.gmt", sep="")
bp.pathways <- gmtPathways(bp.gmt.file)
str(head(bp.pathways))


#### CC
cc.gmt.file <- paste(msigdb_gmt_fp, "c5.go.cc.v2025.1.Hs.entrez.gmt", sep="")
cc.pathways <- gmtPathways(cc.gmt.file)
str(head(cc.pathways))


#### MF

mf.gmt.file <- paste(msigdb_gmt_fp, "c5.go.mf.v2025.1.Hs.entrez.gmt", sep="")
mf.pathways <- gmtPathways(mf.gmt.file)
str(head(mf.pathways))

###########################################################################################
###########################################################################################

# GESECA ALL
setwd("/n/groups/kirschner/Will/BedRest/geseca")

conditions <- colnames(counts)
conditions[grepl("Post", conditions)] <- "Post"
conditions[grepl("Pre", conditions)] <- "Pre"
conditions <- factor(conditions, levels = c("Pre", "Post"))


### Reactome
gr <- geseca(pathways, counts, minSize=15, maxSize=500, scale=TRUE, center=TRUE, eps=0)

jpeg(filename = "reactome.jpg",
     width = 1100,
     height = 1200,
     units = "px",
     pointsize = 12,
     quality = 75 # JPEG specific quality setting
)
plotCoregulationProfile(pathway=pathways[["REACTOME_NEUTROPHIL_DEGRANULATION"]], 
                        E=counts, titles = colnames(counts), conditions=conditions)

dev.off()
plotGesecaTable(gr |> head(10), pathways, E=counts, titles = colnames(counts))

write.csv(gr, "reactome.csv", row.names = FALSE)

###  Hallmark 
gr_hallmark <- geseca(hallmark.pathways, counts, minSize=15, maxSize=500, eps=0)


plotCoregulationProfile(pathway=hallmark.pathways[["HALLMARK_INFLAMMATORY_RESPONSE"]], 
                        E=counts, titles = colnames(counts), conditions=conditions)

plotGesecaTable(gr_hallmark |> head(20), hallmark.pathways, E=counts, titles = colnames(counts))

write.csv(gr_hallmark, "hallmark.csv", row.names = FALSE)

### Wikipathways 
gr_wiki <- geseca(wiki.pathways, counts, minSize=15, maxSize=500, eps=0)


plotCoregulationProfile(pathway=wiki.pathways[["WP_GPCRS_CLASS_A_RHODOPSINLIKE"]], 
                        E=counts, titles = colnames(counts), conditions=conditions)

plotGesecaTable(gr_wiki |> head(20), wiki.pathways, E=counts, titles = colnames(counts))

write.csv(gr_wiki, "wiki.csv", row.names = FALSE)

### KEGG
gr_kegg <- geseca(kegg_pathways, counts, minSize=15, maxSize=500, eps=0)


plotCoregulationProfile(pathway=kegg_pathways[["Cytokine-cytokine receptor interaction"]], 
                        E=counts, titles = colnames(counts), conditions=conditions)

plotGesecaTable(gr_kegg |> head(15), kegg_pathways, E=counts, titles = colnames(counts))
write.csv(gr_kegg, "kegg.csv", row.names = FALSE)

### GO

#### BP
gr_bp <- geseca(bp.pathways, counts, minSize=15, maxSize=500, eps=0)


plotCoregulationProfile(pathway=bp.pathways[["WP_GPCRS_CLASS_A_RHODOPSINLIKE"]], 
                        E=counts, titles = colnames(counts), conditions=conditions)

plotGesecaTable(gr_bp |> head(20), bp.pathways, E=counts, titles = colnames(counts))

write.csv(gr_bp, "bp.csv", row.names = FALSE)
#### CC
gr_cc <- geseca(cc.pathways, counts, minSize=15, maxSize=500, eps=0)


plotCoregulationProfile(pathway=cc.pathways[["WP_GPCRS_CLASS_A_RHODOPSINLIKE"]], 
                        E=counts, titles = colnames(counts), conditions=conditions)

plotGesecaTable(gr_cc |> head(20), cc.pathways, E=counts, titles = colnames(counts))
write.csv(gr_cc, "cc.csv", row.names = FALSE)


#### MF

gr_mf <- geseca(mf.pathways, counts, minSize=15, maxSize=500, eps=0)


plotCoregulationProfile(pathway=mf.pathways[["GOMF_TRANSFERASE_ACTIVITY_TRANSFERRING_ONE_CARBON_GROUPS"]], 
                        E=counts, titles = colnames(counts), conditions=conditions)

plotGesecaTable(gr_mf |> head(20), mf.pathways, E=counts, titles = colnames(counts))

write.csv(gr_mf, "mf.csv", row.names = FALSE)

###########################################################################################
###########################################################################################

# GESECA PRE vs POST
setwd("/n/groups/kirschner/Will/BedRest/geseca")


counts_post <- counts[,grepl("Post", colnames(counts))]
counts_pre <- counts[,grepl("Pre", colnames(counts))]

### Reactome
gr_post <- geseca(pathways, counts_post, minSize=15, maxSize=500)
gr_pre <- geseca(pathways, counts_pre, minSize=15, maxSize=500)
write.csv(gr_post, "reactome_post.csv", row.names = FALSE)
write.csv(gr_pre, "reactome_pre.csv", row.names = FALSE)
###  Hallmark 
gr_post <- geseca(hallmark.pathways, counts_post, minSize=15, maxSize=500)
gr_pre <- geseca(hallmark.pathways, counts_pre, minSize=15, maxSize=500)
write.csv(gr_post, "hallmark_post.csv", row.names = FALSE)
write.csv(gr_pre, "hallmark_pre.csv", row.names = FALSE)

### Wikipathways 
gr_post <- geseca(wiki.pathways, counts_post, minSize=15, maxSize=500)
gr_pre <- geseca(wiki.pathways, counts_pre, minSize=15, maxSize=500)
write.csv(gr_post, "wiki_post.csv", row.names = FALSE)
write.csv(gr_pre, "wiki_pre.csv", row.names = FALSE)

### KEGG
gr_post <- geseca(kegg_pathways, counts_post, minSize=15, maxSize=500)
gr_pre <- geseca(kegg_pathways, counts_pre, minSize=15, maxSize=500)
write.csv(gr_post, "kegg_post.csv", row.names = FALSE)
write.csv(gr_pre, "kegg_pre.csv", row.names = FALSE)

### GO

#### BP
gr_post <- geseca(bp.pathways, counts_post, minSize=15, maxSize=500)
gr_pre <- geseca(bp.pathways, counts_pre, minSize=15, maxSize=500)
write.csv(gr_post, "bp_post.csv", row.names = FALSE)
write.csv(gr_pre, "bp_pre.csv", row.names = FALSE)


#### CC
gr_post <- geseca(cc.pathways, counts_post, minSize=15, maxSize=500)
gr_pre <- geseca(cc.pathways, counts_pre, minSize=15, maxSize=500)
write.csv(gr_post, "cc_post.csv", row.names = FALSE)
write.csv(gr_pre, "cc_pre.csv", row.names = FALSE)



#### MF
gr_post <- geseca(mf.pathways, counts_post, minSize=15, maxSize=500)
gr_pre <- geseca(mf.pathways, counts_pre, minSize=15, maxSize=500)
write.csv(gr_post, "mf_post.csv", row.names = FALSE)
write.csv(gr_pre, "mf_pre.csv", row.names = FALSE)
