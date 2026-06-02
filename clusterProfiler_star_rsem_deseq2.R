library(ggplot2)
library(clusterProfiler)
library(biomaRt)
devtools::load_all("~/projs/enrichplot")
# library(enrichplot)
library(dplyr)
library(stringr)

counts_fp <- "/n/groups/kirschner/Will/BedRest/deseq2/rlog_counts_star_rsem.csv"


counts <- read.csv(counts_fp, sep=",", row.names = 1)
# counts <- counts[,-which(colnames(counts) %in% c("Y11"))]
# counts$transcript_id.s. <- NULL


gene_ids <- sapply(strsplit(rownames(counts), "\\."), head, 1)

# Get the corresponding Entrez IDs
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# attrs <- listAttributes(mart)
genes <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id", "entrezgene_id", "entrezgene_accession"),
               values = gene_ids,
               mart = mart)

nrow(genes)
genes_unique <- genes %>% distinct(ensembl_gene_id, .keep_all = TRUE)
nrow(genes_unique)
genes_unique <- genes_unique %>%
  group_by(entrezgene_id) %>%
  filter(n() == 1) %>%
  ungroup()
nrow(genes_unique)

counts$ensembl_gene_id <- gene_ids
nrow(counts)
counts <- merge(genes_unique, counts, by = "ensembl_gene_id")
nrow(counts)
rownames(counts) <- counts$entrezgene_id
counts$ensembl_gene_id <- NULL
counts$entrezgene_id <- NULL
counts$gene_id <- NULL


############################################
############################################
#### deseq2 results
#
### tmp Post vs Pre
#
deseq2_results_fp <- "/n/groups/kirschner/Will/BedRest/deseq2/star_rsem_tmp_post_pre.csv"
deseq2_results <- read.csv(deseq2_results_fp)
deseq2_results$ensembl_gene_id <- sapply(strsplit(deseq2_results$X, "\\."), head, 1)
deseq2_results <- merge(genes_unique, deseq2_results, by = "ensembl_gene_id")
deseq2_results <- deseq2_results[order(-deseq2_results$stat), ]
# deseq2_results_ <- deseq2_results[!is.na(deseq2_results$stat),]

deseq2_ranks <- deseq2_results$stat
names(deseq2_ranks) <- deseq2_results$entrezgene_id

deseq2_fcs <- deseq2_results$log2FoldChange
names(deseq2_fcs) <- deseq2_results$entrezgene_id

#############################################################################
#############################################################################
##### KEGG
# TODO: https://yulab-smu.top/biomedical-knowledge-mining-book/022-kegg.html#clusterprofiler-kegg-module-gsea

kk <- gseKEGG(geneList     = deseq2_ranks,
              organism     = 'hsa',
              minGSSize    = 15,
              maxGSSize = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE,
              eps=0
              )

# dotplot
dotplot(kk, showCategory=30) #+ ggtitle("dotplot for GSEA")

# cnetplot
edox <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=deseq2_fcs)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=deseq2_fcs)
p3 <- cnetplot(edox, foldChange=deseq2_fcs, circular = TRUE, colorEdge = TRUE) 
plot_list(p1, p2, p3, ncol=3, tag_levels = 'A', widths = c(.8, .8, 1.2))


edox@result$ID_ <- edox@result$ID
edox@result$ID <- edox@result$Description
# treeplot
edox2 <- pairwise_termsim(edox)
edox2@result$setSize <- as.numeric(str_match(edox2@result$leading_edge, "tags=([0-9]+)%")[,2]) / 100

# remove_sets <- c("mitochondrion organization", "")
# xx@compareClusterResult <- filter(xx@compareClusterResult, Description %in% remove_sets)

p <- treeplot(edox2,
              showCategory = 50,
              #size_var = c("Count"), # "setSize"),
              # cluster_panel="dotplot",
              cluster_method = "median",
              color="NES",
              nCluster=10,
              fontsize_tiplab = 3,
              fontsize_cladelab = 4,
              group_color = NULL,
              hexpand=0.1,
              tiplab_offset=6,
              cladelab_offset=70,
              # colnames_angle=45
              # nWords=0,
              # cluster.params = list(method = "ward.D", n = 5, color = "NES", label_words_n = 0, label_format = 30),
              # hilight.params = list(hilight = TRUE, align = "both"),
              # clusterPanel.params = list(clusterPanel = "heatMap", pie = "equal", legend_n = 3,
              # colnames_angle = 0),
              # offset.params = list(bar_tree = 1, cladelab=10, tiplab = 25, extend = 1, hexpand = 2),
) + 
  scale_size_continuous(range = c(1, 4), name="LE %") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0)))

p
treeplot(edox2, cluster_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

# heatplot
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=deseq2_fcs, showCategory=5)
plot_list(p1, p2, ncol=1, tag_levels = 'A')


# output table
kk <- gseKEGG(geneList     = deseq2_ranks,
              organism     = 'hsa',
              minGSSize    = 15,
              maxGSSize    = 1500,
              pvalueCutoff = 1,
              verbose      = TRUE,
              eps=0
)
kk <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')



#############################################################################
#############################################################################
##### Reactome

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ReactomePA")
library(ReactomePA)

react <- gsePathway(deseq2_ranks, 
                    minGSSize = 15,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH", 
                    verbose = TRUE,
                    eps=0
)
head(react)
viewPathway("HCMV Infection", readable=TRUE, foldChange=deseq2_fcs)

# dotplot
dotplot(react, showCategory=30) #+ ggtitle("dotplot for GSEA")

# cnetplot
edox <- setReadable(react, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=deseq2_fcs)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=deseq2_fcs)
p3 <- cnetplot(edox, foldChange=deseq2_fcs, circular = TRUE, colorEdge = TRUE) 
plot_list(p1, p2, p3, ncol=3, tag_levels = 'A', widths = c(.8, .8, 1.2))

# treeplot
edox@result$ID_ <- edox@result$ID
edox@result$ID <- edox@result$Description
# treeplot
edox2 <- pairwise_termsim(edox)
edox2@result$setSize <- as.numeric(str_match(edox2@result$leading_edge, "tags=([0-9]+)%")[,2]) / 100

# remove_sets <- c("mitochondrion organization", "")
# xx@compareClusterResult <- filter(xx@compareClusterResult, Description %in% remove_sets)

p <- treeplot(edox2,
              showCategory = 50,
              #size_var = c("Count"), # "setSize"),
              # cluster_panel="dotplot",
              cluster_method = "average",
              color="NES",
              nCluster=10,
              fontsize_tiplab = 3,
              fontsize_cladelab = 4,
              group_color = NULL,
              hexpand=0.1,
              tiplab_offset=6,
              cladelab_offset=70,
              # colnames_angle=45
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

p

# heatplot
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=deseq2_fcs, showCategory=5)
plot_list(p1, p2, ncol=1, tag_levels = 'A')


# output table
library(ReactomePA)
react <- gsePathway(deseq2_ranks, 
                    minGSSize = 15,
                    maxGSSize = 1500,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH", 
                    verbose = TRUE,
                    eps=0
)
react <- setReadable(react, 'org.Hs.eg.db', 'ENTREZID')

react_age <- gsePathway(geneList     = deseq2_ranks_age,
                        minGSSize    = 15,
                        maxGSSize    = 1500,
                        pvalueCutoff = 1,
                        verbose      = TRUE,
                        pAdjustMethod = "BH", 
                        eps=0
)
react_age <- setReadable(react_age, 'org.Hs.eg.db', 'ENTREZID')

react_merged <- full_join(react@result, react_age@result, by = "ID", suffix = c("_virus","_age"))
nrow(react_merged)
react_merged <- react_merged[((!is.na(react_merged$p.adjust_virus) & (react_merged$p.adjust_virus < 0.05))) | (!is.na(react_merged$p.adjust_age) & (react_merged$p.adjust_age < 0.05)),]

react_merged$Description <- react_merged$Description_age
react_merged$Description[is.na(react_merged$Description)] <- react_merged$Description_virus[is.na(react_merged$Description)]
react_merged$setSize <- react_merged$setSize_age
react_merged$setSize[is.na(react_merged$setSize)] <- react_merged$setSize_virus[is.na(react_merged$setSize)]

react_merged$core_enrichment_virus <- str_replace_all(react_merged$core_enrichment_virus, "/", ", ")
react_merged$core_enrichment_age <- str_replace_all(react_merged$core_enrichment_age, "/", ", ")
react_merged <- react_merged[c("ID", "Description", "setSize", "NES_virus", "p.adjust_virus", "core_enrichment_virus", "NES_age", "p.adjust_age", "core_enrichment_age")]

# react_merged$NES_virus <- replace_na_0(react_merged$NES_virus)
# react_merged$NES_age <- replace_na_0(react_merged$NES_age)
# react_merged$p.adjust_virus <- replace_na_1(react_merged$p.adjust_virus)
# react_merged$p.adjust_age <- replace_na_1(react_merged$p.adjust_age)

setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2/fgsea")
write.table(react_merged, "star_rsem_clusterProfiler_virus_age_Reactome.csv", row.names = FALSE, sep="\t", quote=FALSE)




#############################################################################
#############################################################################
##### GO BP
library(org.Hs.eg.db)

go_bp_levels <- read.csv("/n/groups/kirschner/reference/GO/go_bp_id_to_name.csv")
go_bp_levels <- go_bp_levels %>% distinct(id, name, .keep_all = TRUE)
rownames(go_bp_levels) <- go_bp_levels$id


ego <- gseGO(geneList     = deseq2_ranks,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP",
             minGSSize    = 15,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = TRUE,
             eps=0
)
# ego@result <- ego@result %>%
#   filter(ID %in% go_bp_levels_5$id)

nrow(ego@result)

# dotplot
p1 <- dotplot(ego, showCategory=30) #+ ggtitle("dotplot for ego")
p2 <- dotplot(ego1, showCategory=30) + ggtitle("dotplot for ego1")
plot_list(p1, p2)
p1
# cnetplot
ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(ego, foldChange=deseq2_fcs)
p1
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(ego, categorySize="pvalue", foldChange=deseq2_fcs)
p2
p3 <- cnetplot(edox, foldChange=deseq2_fcs, circular = TRUE, colorEdge = TRUE) 
plot_list(p1, p2, p3, ncol=3, tag_levels = 'A', widths = c(.8, .8, 1.2))

ego <- arrange(ego, p.adjust)
# padj <- ego$p.adjust
# padj[ego$NES<0] <- -padj[ego$NES<0]
# ego$p.adjust <- padj


# emapplot
egox <- pairwise_termsim(ego)
p1 <- emapplot(egox, color="NES")
p1

egox2 <-simplify(egox, cutoff=0.5, by="p.adjust", select_fun=min)
p2 <- emapplot(egox2)
egox3 <-simplify(egox, cutoff=0.9, by="p.adjust", select_fun=min)
p3 <- emapplot(egox3)

plot_list(p1, p2, p3, ncol=3, widths = c(1.2, 1, .8))
plot_list(p1, p3, ncol=2)
p3

ego_pos <- filter(ego, p.adjust < 0.05 & NES > 0)
ego_neg <- filter(ego, p.adjust < 0.05 & NES < 0)

# treeplot
egox <- pairwise_termsim(ego)
p1 <- treeplot(egox, showCategory = 70, color="NES", extend=0.1) #, cluster.params = list(n = 10))
p1
p1 + scale_color_gradient(low = "deepskyblue2", high = "brown2", name="NES")
p1
# treeplot(egox, showCategory = 150, color="NES")

ego2 <-simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
egox2 <- pairwise_termsim(ego2)
p2 <- treeplot(egox2, showCategory = 50, color="NES", cluster.params = list(n = 5))
p2 + scale_color_gradient(low = "deepskyblue2", high = "brown2", name="NES")

ego3 <-simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
egox3 <- pairwise_termsim(ego3)
p3 <- treeplot(egox3, showCategory = 50, color="NES", cluster.params = list(n = 6))
p3

ego4 <-simplify(ego, cutoff=0.5, by="p.adjust", select_fun=min)
egox4 <- pairwise_termsim(ego4)
p4 <- treeplot(egox4, showCategory = 50, color="NES")

plot_list(p1, p2, p3, p4, ncol=2)
plot_list(p1, p3, ncol=2)
plot_list(p3, p4, ncol=2)
p3

treeplot(edox2, showCategory = 150)
treeplot(edox2, cluster_method = "average", showCategory = 50)
# aplot::plot_list(p1, p2, tag_levels='A')

########################

go_bp_levels <- read.csv("/n/groups/kirschner/reference/GO/go_bp_id_to_name.csv")
go_bp_levels <- go_bp_levels %>% distinct(id, name, .keep_all = TRUE)
rownames(go_bp_levels) <- go_bp_levels$id

ego@result$level <- go_bp_levels[ego@result$ID, "level"]
# setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2")
# write.csv(xx@compareClusterResult, "age_vs_virus_GO_BP.csv")

ego@result <- filter(ego@result, level > 3)
# ego@result[xx@compareClusterResult$ID == "GO:0002460", "Description"] <- "adaptive immune response based on somatic recombination of immune receptors*"

ego@result$ID_ <- ego@result$ID
ego@result$ID <- ego@result$Description

xx <- pairwise_termsim(ego)       


xx@result$setSize <- as.numeric(str_match(xx@result$leading_edge, "tags=([0-9]+)%")[,2]) / 100

# remove_sets <- c("mitochondrion organization", "")
# xx@compareClusterResult <- filter(xx@compareClusterResult, Description %in% remove_sets)

treeplot(xx,
         showCategory = 80,
         #size_var = c("Count"), # "setSize"),
         # cluster_panel="dotplot",
         cluster_method = "median",
         color="NES",
         nCluster=7,
         fontsize_tiplab = 4,
         fontsize_cladelab = 4,
         group_color = NULL,
         hexpand=0.1,
         tiplab_offset=4,
         cladelab_offset=70,
         # colnames_angle=45
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

xx2 <- clusterProfiler::simplify(xx, cutoff=0.7, by="setSize", select_fun=max)

xx2 <- pairwise_termsim(xx2)

p <- treeplot(xx2,
              showCategory = 80,
              #size_var = c("Count"), # "setSize"),
              cluster_panel="dotplot",
              cluster_method = "complete",
              color="NES",
              nCluster=7,
              fontsize_tiplab = 4,
              fontsize_cladelab = 4,
              group_color = NULL,
              hexpand=0.05,
              tiplab_offset=4,
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

p
setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2")
ggsave("age_vs_virus.png", plot = p, dpi = 400)

# dotplot(xx, color="NES", showCategory = 35)
# use `geneClusterPanel` to change the gene cluster panel.

p <- treeplot(
  xx,
  geneClusterPanel = "dotplot",
  size = "setSize",
  cluster.params = list(method = "ward.D", n = 7, label_format = 30),
  showCategory = 70,
  color = "NES",
  offset.params = list(bar_tree=rel(2), tiplab = 4, hexpand = 0.5)
) + 
  scale_color_gradient(low = "deepskyblue2", high = "brown2", name="NES") #+
# scale_size_continuous(range = c(1.8, 5), name="size")
p$data$Count <- p$data$count
p
xx@compareClusterResult$setSize <- as.numeric(str_match(xx@compareClusterResult$leading_edge, "tags=([0-9]+)%")[,2]) / 100

xx2 <-simplify(xx, cutoff=0.9, by="setSize", select_fun=max)
xx2 <- pairwise_termsim(xx2)  
treeplot(
  xx2,
  geneClusterPanel = "dotplot",
  size = "setSize",
  cluster.params = list(method = "ward.D", n = 7, label_format = 30),
  showCategory = 60,
  color = "NES",
  offset.params = list(bar_tree=rel(2), tiplab = 4, hexpand = 0.5)
) + 
  scale_color_gradient(low = "deepskyblue2", high = "brown2", name="NES")

# treeplot(xx, geneClusterPanel = "pie")  
p <- goplot(ego)
p
#############################################################################
#############################################################################
##### GO CC

ego_cc <- gseGO(geneList     = deseq2_ranks,
                OrgDb        = org.Hs.eg.db,
                ont          = "CC",
                minGSSize    = 15,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE,
                eps=0
                )
head(ego_cc)

# dotplot
dotplot(ego_cc, showCategory=30) #+ ggtitle("dotplot for GSEA")

# cnetplot
edox <- setReadable(ego_cc, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=deseq2_fcs)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=deseq2_fcs)
p3 <- cnetplot(edox, foldChange=deseq2_fcs, circular = TRUE, colorEdge = TRUE) 
plot_list(p1, p2, p3, ncol=3, tag_levels = 'A', widths = c(.8, .8, 1.2))

# emapplot
egox <- pairwise_termsim(egox)
p1 <- emapplot(egox)
egox2 <-simplify(egox, cutoff=0.5, by="p.adjust", select_fun=min)
p2 <- emapplot(egox2)
egox3 <-simplify(egox, cutoff=0.9, by="p.adjust", select_fun=min)
p3 <- emapplot(egox3)

plot_list(p1, p2, p3, ncol=3, widths = c(1.2, 1, .8))
plot_list(p1, p3, ncol=2)
p3

# treeplot
p1 <- treeplot(egox)
treeplot(edox2, cluster_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

# heatplot
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=deseq2_fcs, showCategory=5)
plot_list(p1, p2, ncol=1, tag_levels = 'A')

library(org.Hs.eg.db)

go_cc_levels <- read.csv("/n/groups/kirschner/reference/GO/go_cc_id_to_name.csv")
go_cc_levels <- go_cc_levels %>% distinct(id, name, .keep_all = TRUE)
rownames(go_cc_levels) <- go_cc_levels$id

ego <- gseGO(geneList     = deseq2_ranks,
             OrgDb        = org.Hs.eg.db,
             ont          = "CC",
             minGSSize    = 15,
             maxGSSize    = 1500,
             pvalueCutoff = 1,
             verbose      = TRUE,
             eps=0
)
ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')

ego_age <- gseGO(geneList     = deseq2_ranks_age,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "CC",
                 minGSSize    = 15,
                 maxGSSize    = 1500,
                 pvalueCutoff = 1,
                 verbose      = TRUE,
                 eps=0
)
ego_age <- setReadable(ego_age, 'org.Hs.eg.db', 'ENTREZID')

ego_merged <- full_join(ego@result, ego_age@result, by = "ID", suffix = c("_virus","_age"))
nrow(ego_merged)
ego_merged <- ego_merged[((!is.na(ego_merged$p.adjust_virus) & (ego_merged$p.adjust_virus < 0.05))) | (!is.na(ego_merged$p.adjust_age) & (ego_merged$p.adjust_age < 0.05)),]

nrow(ego@result)
nrow(ego_age@result)
nrow(ego_merged)
ego_merged$level <- go_cc_levels[ego_merged$ID, "level"]
ego_merged$msigdb_name <- go_cc_levels[ego_merged$ID, "msigdb_name"]

ego_merged$Description <- ego_merged$Description_age
ego_merged$Description[is.na(ego_merged$Description)] <- ego_merged$Description_virus[is.na(ego_merged$Description)]
ego_merged$setSize <- ego_merged$setSize_age
ego_merged$setSize[is.na(ego_merged$setSize)] <- ego_merged$setSize_virus[is.na(ego_merged$setSize)]

ego_merged$core_enrichment_virus <- str_replace_all(ego_merged$core_enrichment_virus, "/", ", ")
ego_merged$core_enrichment_age <- str_replace_all(ego_merged$core_enrichment_age, "/", ", ")
ego_merged <- ego_merged[c("ID", "level", "Description", "setSize", "NES_virus", "p.adjust_virus", "core_enrichment_virus", "NES_age", "p.adjust_age", "core_enrichment_age", "msigdb_name")]

# ego_merged$NES_virus <- replace_na_0(ego_merged$NES_virus)
# ego_merged$NES_age <- replace_na_0(ego_merged$NES_age)
# ego_merged$p.adjust_virus <- replace_na_1(ego_merged$p.adjust_virus)
# ego_merged$p.adjust_age <- replace_na_1(ego_merged$p.adjust_age)

setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2/fgsea")
write.table(ego_merged, "star_rsem_clusterProfiler_virus_age_GOCC.csv", row.names = FALSE, sep="\t", quote=FALSE)

#############################################################################
#############################################################################
##### GO MF

ego_mf <- gseGO(geneList     = deseq2_ranks,
                OrgDb        = org.Hs.eg.db,
                ont          = "MF",
                minGSSize    = 15,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
head(ego_mf)

# dotplot
dotplot(ego_mf, showCategory=30) #+ ggtitle("dotplot for GSEA")

# cnetplot
edox <- setReadable(ego_mf, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=deseq2_fcs)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=deseq2_fcs)
p3 <- cnetplot(edox, foldChange=deseq2_fcs, circular = TRUE, colorEdge = TRUE) 
plot_list(p1, p2, p3, ncol=3, tag_levels = 'A', widths = c(.8, .8, 1.2))

# treeplot
edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
treeplot(edox2, cluster_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

# heatplot
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=deseq2_fcs, showCategory=5)
plot_list(p1, p2, ncol=1, tag_levels = 'A')


library(org.Hs.eg.db)

go_mf_levels <- read.csv("/n/groups/kirschner/reference/GO/go_mf_id_to_name.csv")
go_mf_levels <- go_mf_levels %>% distinct(id, name, .keep_all = TRUE)
rownames(go_mf_levels) <- go_mf_levels$id

ego <- gseGO(geneList     = deseq2_ranks,
             OrgDb        = org.Hs.eg.db,
             ont          = "MF",
             minGSSize    = 15,
             maxGSSize    = 1500,
             pvalueCutoff = 1,
             verbose      = TRUE,
             eps=0
)
ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')

ego_age <- gseGO(geneList     = deseq2_ranks_age,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "MF",
                 minGSSize    = 15,
                 maxGSSize    = 1500,
                 pvalueCutoff = 1,
                 verbose      = TRUE,
                 eps=0
)
ego_age <- setReadable(ego_age, 'org.Hs.eg.db', 'ENTREZID')

ego_merged <- full_join(ego@result, ego_age@result, by = "ID", suffix = c("_virus","_age"))
nrow(ego_merged)
ego_merged <- ego_merged[((!is.na(ego_merged$p.adjust_virus) & (ego_merged$p.adjust_virus < 0.05))) | (!is.na(ego_merged$p.adjust_age) & (ego_merged$p.adjust_age < 0.05)),]

nrow(ego@result)
nrow(ego_age@result)
nrow(ego_merged)
ego_merged$level <- go_mf_levels[ego_merged$ID, "level"]
ego_merged$msigdb_name <- go_mf_levels[ego_merged$ID, "msigdb_name"]

ego_merged$Description <- ego_merged$Description_age
ego_merged$Description[is.na(ego_merged$Description)] <- ego_merged$Description_virus[is.na(ego_merged$Description)]
ego_merged$setSize <- ego_merged$setSize_age
ego_merged$setSize[is.na(ego_merged$setSize)] <- ego_merged$setSize_virus[is.na(ego_merged$setSize)]

ego_merged$core_enrichment_virus <- str_replace_all(ego_merged$core_enrichment_virus, "/", ", ")
ego_merged$core_enrichment_age <- str_replace_all(ego_merged$core_enrichment_age, "/", ", ")
ego_merged <- ego_merged[c("ID", "level", "Description", "setSize", "NES_virus", "p.adjust_virus", "core_enrichment_virus", "NES_age", "p.adjust_age", "core_enrichment_age", "msigdb_name")]

# ego_merged$NES_virus <- replace_na_0(ego_merged$NES_virus)
# ego_merged$NES_age <- replace_na_0(ego_merged$NES_age)
# ego_merged$p.adjust_virus <- replace_na_1(ego_merged$p.adjust_virus)
# ego_merged$p.adjust_age <- replace_na_1(ego_merged$p.adjust_age)

setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2/fgsea")
write.table(ego_merged, "star_rsem_clusterProfiler_virus_age_GOMF.csv", row.names = FALSE, sep="\t", quote=FALSE)

#############################################################################
#############################################################################
##### Wiki

wiki <- gseWP(deseq2_ranks,
              organism = "Homo sapiens",
              minGSSize    = 15,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
head(wiki)

# dotplot
dotplot(wiki, showCategory=30) + ggtitle("dotplot for GSEA")

# cnetplot
edox <- setReadable(wiki, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=deseq2_fcs)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=deseq2_fcs)
p3 <- cnetplot(edox, foldChange=deseq2_fcs, circular = TRUE, colorEdge = TRUE) 
plot_list(p1, p2, p3, ncol=3, tag_levels = 'A', widths = c(.8, .8, 1.2))

# treeplot
edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
treeplot(edox2, cluster_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

# heatplot
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=deseq2_fcs, showCategory=5)
plot_list(p1, p2, ncol=1, tag_levels = 'A')

# output table
wiki <- gseWP(deseq2_ranks,
              organism = "Homo sapiens",
              minGSSize    = 15,
              maxGSSize    = 1500,
              pvalueCutoff = 1,
              verbose      = TRUE,
              eps=0
)

wiki <- setReadable(wiki, 'org.Hs.eg.db', 'ENTREZID')

wiki_age <- gseWP(deseq2_ranks_age,
                  organism = "Homo sapiens",
                  minGSSize    = 15,
                  maxGSSize    = 1500,
                  pvalueCutoff = 1,
                  verbose      = TRUE,
                  eps=0
)

wiki_age <- setReadable(wiki_age, 'org.Hs.eg.db', 'ENTREZID')

wiki_merged <- full_join(wiki@result, wiki_age@result, by = "ID", suffix = c("_virus","_age"))
nrow(wiki_merged)
wiki_merged <- wiki_merged[((!is.na(wiki_merged$p.adjust_virus) & (wiki_merged$p.adjust_virus < 0.05))) | (!is.na(wiki_merged$p.adjust_age) & (wiki_merged$p.adjust_age < 0.05)),]

wiki_merged$Description <- wiki_merged$Description_age
wiki_merged$Description[is.na(wiki_merged$Description)] <- wiki_merged$Description_virus[is.na(wiki_merged$Description)]
wiki_merged$setSize <- wiki_merged$setSize_age
wiki_merged$setSize[is.na(wiki_merged$setSize)] <- wiki_merged$setSize_virus[is.na(wiki_merged$setSize)]

wiki_merged$core_enrichment_virus <- str_replace_all(wiki_merged$core_enrichment_virus, "/", ", ")
wiki_merged$core_enrichment_age <- str_replace_all(wiki_merged$core_enrichment_age, "/", ", ")
wiki_merged <- wiki_merged[c("ID", "Description", "setSize", "NES_virus", "p.adjust_virus", "core_enrichment_virus", "NES_age", "p.adjust_age", "core_enrichment_age")]

# wiki_merged$NES_virus <- replace_na_0(wiki_merged$NES_virus)
# wiki_merged$NES_age <- replace_na_0(wiki_merged$NES_age)
# wiki_merged$p.adjust_virus <- replace_na_1(wiki_merged$p.adjust_virus)
# wiki_merged$p.adjust_age <- replace_na_1(wiki_merged$p.adjust_age)

setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2/fgsea")
write.table(wiki_merged, "star_rsem_clusterProfiler_virus_age_WikiPathways.csv", row.names = FALSE, sep="\t", quote=FALSE)



#############################################################################
#############################################################################
#### Hallmark
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens")

H_t2g <- m_df |>
  dplyr::filter(gs_collection == "H") |>
  dplyr::select(gs_name, ncbi_gene)

hm <- GSEA(deseq2_ranks,
           TERM2GENE = H_t2g,
           minGSSize    = 15,
           maxGSSize    = 1500,
           pvalueCutoff = 1,
           verbose      = TRUE,
           eps=0
)

hm <- setReadable(hm, 'org.Hs.eg.db', 'ENTREZID')

hm_age <- GSEA(deseq2_ranks_age,
               TERM2GENE = H_t2g,
               minGSSize    = 15,
               maxGSSize    = 1500,
               pvalueCutoff = 1,
               verbose      = TRUE,
               eps=0
)

hm_age <- setReadable(hm_age, 'org.Hs.eg.db', 'ENTREZID')

hm_merged <- full_join(hm@result, hm_age@result, by = "ID", suffix = c("_virus","_age"))
nrow(hm_merged)
hm_merged <- hm_merged[((!is.na(hm_merged$p.adjust_virus) & (hm_merged$p.adjust_virus < 0.05))) | (!is.na(hm_merged$p.adjust_age) & (hm_merged$p.adjust_age < 0.05)),]

hm_merged$Description <- hm_merged$Description_age
hm_merged$Description[is.na(hm_merged$Description)] <- hm_merged$Description_virus[is.na(hm_merged$Description)]
hm_merged$setSize <- hm_merged$setSize_age
hm_merged$setSize[is.na(hm_merged$setSize)] <- hm_merged$setSize_virus[is.na(hm_merged$setSize)]

hm_merged$core_enrichment_virus <- str_replace_all(hm_merged$core_enrichment_virus, "/", ", ")
hm_merged$core_enrichment_age <- str_replace_all(hm_merged$core_enrichment_age, "/", ", ")
hm_merged <- hm_merged[c("ID", "Description", "setSize", "NES_virus", "p.adjust_virus", "core_enrichment_virus", "NES_age", "p.adjust_age", "core_enrichment_age")]

# hm_merged$NES_virus <- replace_na_0(hm_merged$NES_virus)
# hm_merged$NES_age <- replace_na_0(hm_merged$NES_age)
# hm_merged$p.adjust_virus <- replace_na_1(hm_merged$p.adjust_virus)
# hm_merged$p.adjust_age <- replace_na_1(hm_merged$p.adjust_age)

setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2/fgsea")
write.table(hm_merged, "star_rsem_clusterProfiler_virus_age_Hallmark.csv", row.names = FALSE, sep="\t", quote=FALSE)





################################################################################
################################################################################
################################################################################
# Treeplot 
##########
library(GO.db)
library(org.Hs.eg.db)

setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2/fgsea")
gobp_filter <- read.csv("GSEA_adipose_virus_and_age_GOBP.csv")
gobp_filter[gobp_filter$Description == "neuron differentiation", "comment"] = 0
gobp_filter[gobp_filter$Description == "cell morphogenesis involved in neuron differentiation", "comment"] = 0
gobp_filter[gobp_filter$Description == "tube development", "comment"] = 0
gobp_filter[gobp_filter$Description == "developmental growth involved in morphogenesis", "comment"] = 0
gobp_filter[gobp_filter$Description == "cell junction organization", "comment"] = 0

kegg_filter <- read.csv("GSEA_adipose_virus_and_age_KEGG.csv")
reactome_filter <- read.csv("GSEA_adipose_virus_and_age_Reactome.csv")

#### combined analysis
## KEGG
kegg.data <- download_KEGG(species="hsa", keggType = "KEGG", keyType = "kegg")
term2gene.kegg <- kegg.data$KEGGPATHID2EXTID
term2name.kegg <- kegg.data$KEGGPATHID2NAME
term2name.kegg$to <- tolower(term2name.kegg$to)
term2name.kegg$to <- paste0(term2name.kegg$to, " (KEGG)")

kegg_filter_leave <- kegg_filter %>% filter(comment == 1)
term2gene.kegg <- term2gene.kegg %>% filter(from %in% kegg_filter_leave$ID)

## GO BP

ont <- "BP"  ##limit sets to GO-BP
goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
if (ont != "ALL") {goterms <- goterms[goterms == ont]}
term2gene.go <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=names(goterms), column="ENTREZID", keytype="GOALL", multiVals='list')
term2gene.go <- stack(term2gene.go)
term2gene.go <- term2gene.go[!is.na(term2gene.go[,"values"]),] [,c(2,1)]
colnames(term2gene.go) <- c("from","to")

term2name.go <- Term(GOTERM)[ names(goterms) ]
term2name.go <- data.frame("from"=names(term2name.go), "to"=term2name.go)

setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2/fgsea")
write.csv(term2gene.go, "term2gene.go.csv", )
write.csv(term2name.go, "term2name.go.csv")

gobp_filter_remove <- gobp_filter %>% filter(comment == 0)
nrow(term2gene.go)
term2gene.go <- term2gene.go %>% filter(!(from %in% gobp_filter_remove$ID))
nrow(term2gene.go)

go_bp_levels <- read.csv("/n/groups/kirschner/reference/GO/go_bp_id_to_name.csv")
go_bp_levels <- go_bp_levels %>% distinct(id, name, .keep_all = TRUE)
rownames(go_bp_levels) <- go_bp_levels$id

go_bp_levels_456 <- go_bp_levels[(go_bp_levels$level >= 4) & (go_bp_levels$level <= 6),]

ego <- gseGO(geneList     = deseq2_ranks,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP",
             minGSSize    = 15,
             maxGSSize    = 1500,
             pvalueCutoff = 0.05,
             verbose      = TRUE,
             eps=0
)
ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')

ego@result <- ego@result %>%
  filter(ID %in% go_bp_levels_456$id)

ego@result <- ego@result %>%
  filter(!(ID %in% gobp_filter_remove$ID))

ego@result$setSize <- as.numeric(str_match(ego@result$leading_edge, "tags=([0-9]+)%")[,2]) / 100


ego <- arrange(ego, p.adjust)

ego@result$level <- go_bp_levels[ego@result$ID,]

ego <- pairwise_termsim(ego) #, showCategory = 200)
p1 <- treeplot(ego,
               showCategory = 70,
               color = "NES",
               # size_var = c("Count", "setSize"),
               nCluster = 6,
               cluster_method = "ward.D",
               label_format = 30,
               fontsize_tiplab = 4,
               fontsize_cladelab = 4,
               group_color = NULL,
               extend = 0.3,
               hilight = TRUE,
               align = "both",
               hexpand = 0.1,
               tiplab_offset = 0.2,
               cladelab_offset = 10
)
p1 + scale_color_gradient(low = "deepskyblue2", high = "brown2", name="NES")+
  scale_size_continuous(range = c(1.8, 5), name="LE %")

tree_data <- as.data.frame(p1@data)

# ego2 <- simplify(ego, cutoff=0.67, by="p.adjust", select_fun=min)

# View(ego2@result)

ego@result$setSize <- as.numeric(str_match(ego2@result$leading_edge, "tags=([0-9]+)%")[,2]) / 100

ego2 <- simplify(ego, cutoff=0.8, by="setSize", select_fun=max)

ego2 <- pairwise_termsim(ego2) #, showCategory = 200)

ego2@result[ego2@result$ID == "GO:0002460", "Description"] <- "adaptive immune response based on somatic recombination of immune receptors*"
# built from immunoglobulin superfamily domains

library(ggtree)

p2 <- treeplot(ego2,
               showCategory = 60,
               color = "NES",
               # size_var = c("leading_edge_tags"),
               nCluster = 7,
               cluster_method = "ward.D",
               label_format = 30,#function(x) {str_wrap(x, width = 20)},
               fontsize_tiplab = 3,
               fontsize_cladelab = 4,
               group_color = NULL,
               extend = 0.3,
               hilight = TRUE,
               align = "both",
               hexpand = 0.2,
               tiplab_offset = 0.3,
               cladelab_offset = 20
) + scale_color_gradient(low = "deepskyblue2", high = "brown2", name="NES") +
  scale_size_continuous(range = c(1.8, 5), name="LE %") #+
# geom_tiplab(aes(label=str_wrap(label, width=70)), offset=0.3, size=3, geom="text") #, label.padding=0)

p2
setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2")
ggsave("treeplot.png", plot = p2, dpi = 400)

library(tidytree)
tree_object <- as.treedata(p2@data)

## Reactome
library(fgsea)
msigdb_gmt_fp <- "/n/groups/kirschner/reference/msigdb/msigdb_v2025.1.Hs_files_to_download_locally/msigdb_v2025.1.Hs_GMTs/"
gmt.file <- paste(msigdb_gmt_fp, "c2.cp.reactome.v2025.1.Hs.entrez.gmt", sep="")
pathways <- gmtPathways(gmt.file)
str(head(pathways))

library(purrr)

term2gene.reactome <- imap_dfr(pathways, function(x, y) {data.frame("from"=y, "to"=x)})
reactome_names <-  tolower(gsub("_", " " ,sub("REACTOME_", "", names(pathways))))
term2name.reactome <- data.frame("from"=names(pathways), "to"=reactome_names)
term2name.reactome$to <- paste0(term2name.reactome$to, " (Reactome)")

reactome_filter$NAME <- paste0("REACTOME_", gsub(" ", "_", toupper(reactome_filter$Description)))
reactome_filter_leave <- reactome_filter %>% filter(comment == 1)
reactome_filter_leave <- reactome_filter_leave %>% add_row(NAME="REACTOME_VIRAL_MESSENGER_RNA_SYNTHESIS")
term2gene.reactome <- term2gene.reactome %>% filter(from %in% reactome_filter_leave$NAME)


## combine the KEGG and GO TERM2GENE and TERM2NAME
TERM2GENE <- rbind(term2gene.kegg, term2gene.reactome)
TERM2NAME <- rbind(term2name.kegg, term2name.reactome)

res.combined <- GSEA(
  geneList = deseq2_ranks,
  minGSSize = 10,
  maxGSSize = 2500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  eps = 0
)
res.combined <- setReadable(res.combined, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

res.combined <- pairwise_termsim(res.combined)       

# res.combined2 <- clusterProfiler::simplify(res.combined, cutoff=0.7, by="p.adjust", select_fun=min)

treeplot(res.combined,
         showCategory = 50,
         color = "p.adjust",
         # size_var = c("Count", "setSize"),
         nCluster = 2,
         cluster_method = "ward.D",
         label_format = 30,
         fontsize_tiplab = 4,
         fontsize_cladelab = 0,
         group_color = c("black", "black"),
         extend = 0.3,
         hilight = TRUE,
         align = "both",
         hexpand = 0.1,
         tiplab_offset = 0.3,
         cladelab_offset = 14
) + scale_size_continuous(range = c(2, 6), name="setSize") + scale_fill_manual(values = c("white", "white"))
# label_format 	wrap length for labels or custom formatting function
# fontsize_tiplab 	font size for tip labels
# fontsize_cladelab 	font size for clade labels
# group_color 	vector of colors for groups
# extend 	extend length for clade labels
# hilight 	whether to highlight clades
# align 	alignment for highlight rectangles
# hexpand 	expand x limits by amount of xrange * hexpand
# tiplab_offset 	offset for tip labels
# cladelab_offset offset for clade labels


virus_kegg <- term2name.kegg[(grepl("virus", term2name.kegg$to)) & (grepl("infection", term2name.kegg$to)),]
virus_kegg <- virus_kegg[virus_kegg$to != "human papillomavirus infection (KEGG)",]
virus_kegg <- virus_kegg[virus_kegg$to != "human t-cell leukemia virus 1 infection (KEGG)",]
virus_kegg_genes <- term2gene.kegg[term2gene.kegg$from %in% virus_kegg$from,]

res.virus <- GSEA(
  geneList = deseq2_ranks,
  minGSSize = 1,
  maxGSSize = 2500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = virus_kegg_genes,
  TERM2NAME = TERM2NAME,
  eps = 0
)

virus_gene_sets <- split(virus_kegg_genes$to, virus_kegg_genes$from)
common_genes <- Reduce(intersect, virus_gene_sets)


genes_unique_ <- as.data.frame(genes_unique)
rownames(genes_unique_) <- genes_unique_$entrezgene_id

common_genes_names <- genes_unique_[common_genes, "entrezgene_accession"]

# virus_kegg_genes %>%
#   group_by(from) %>%
#   summarise(values = list(to)) %>%
#   pull(values)

################################################################################
################################################################################
################################################################################
### Heatmaps

# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

samples_fp <- '/n/groups/kirschner/Will/HCMV/CMV_body_characteristics.csv'
samples <- read.csv(samples_fp, sep=",")

rownames(samples) <- samples$participant_id
samples <- samples[which(rownames(samples) %in% colnames(counts)),]
samples$age <- factor(samples$age_group)
samples$virus <- "CMV-"
samples[samples$CMV != "<< Y range", "virus"] <- "CMV+"
samples$virus <- factor(samples$virus, levels=c("CMV-", "CMV+"))
samples$age <- factor(samples$age)
levels(samples$age) <- c("Old", "Young")

samples <- samples[order(samples$virus, samples$age), ]

samples[samples$CMV == "<< Y range", "CMV"] <- 0
samples$CMV <- as.numeric(samples$CMV)


row_ha = rowAnnotation(
  age = samples$age,
  virus = samples$virus,
  "CMV IgG (IU/mL)" = anno_barplot(samples$CMV, axis_gp = gpar(fontsize = 15)),
  col=list(virus=c("CMV-" = "#d4d4d4", "CMV+" = "#b256f5"), age=c("Yong"="#2dcc30", "Old"="#f5ef38")),
  gap = unit(c(1, 2), "mm")
)

ht_opt$message = FALSE

genes_unique <- as.data.frame(genes_unique)
rownames(genes_unique) <- genes_unique$entrezgene_id

################################################################################
### GO BP
library(org.Hs.eg.db)

ego <- gseGO(geneList     = deseq2_ranks,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP",
             minGSSize    = 15,
             maxGSSize    = 1500,
             pvalueCutoff = 1,
             verbose      = TRUE,
             eps=0
)
ego_age <- gseGO(geneList     = deseq2_ranks_age,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 15,
                 maxGSSize    = 1500,
                 pvalueCutoff = 1,
                 verbose      = TRUE,
                 eps=0
)

ego_merged <- full_join(ego@result, ego_age@result, by = "ID", suffix = c("_virus","_age"))
ego_merged <- ego_merged[((!is.na(ego_merged$p.adjust_virus) & (ego_merged$p.adjust_virus < 0.05))) | (!is.na(ego_merged$p.adjust_age) & (ego_merged$p.adjust_age < 0.05)),]
ego_merged$Description <- ego_merged$Description_age
ego_merged$Description[is.na(ego_merged$Description)] <- ego_merged$Description_virus[is.na(ego_merged$Description)]
ego_merged$setSize <- ego_merged$setSize_age
ego_merged$setSize[is.na(ego_merged$setSize)] <- ego_merged$setSize_virus[is.na(ego_merged$setSize)]

ego_merged$core_enrichment_virus <- strsplit(ego_merged$core_enrichment_virus, "/")
ego_merged$core_enrichment_age <- strsplit(ego_merged$core_enrichment_age, "/")

pathway <- "response to virus"
n_top = 70

col <- "core_enrichment_virus"

pathway_genes <- ego_merged[ego_merged$Description == pathway,][col][[1]][[1]][1:n_top]
pathway_genes <- pathway_genes[pathway_genes %in% rownames(counts)]
gene_counts <- counts[pathway_genes, samples$participant_id]
gene_counts <- t(as.matrix(gene_counts))
gene_counts <- scale(gene_counts)
gene_counts <- gene_counts[,!apply(gene_counts, 2, function(x) any(is.nan(x)))]
colnames(gene_counts) <- genes_unique[colnames(gene_counts), "entrezgene_accession"]

row_ha = rowAnnotation(
  age = samples$age,
  virus = samples$virus,
  "CMV IgG\n(IU/mL)" = anno_barplot(samples$CMV),
  col=list(virus=c("CMV-" = "#d4d4d4", "CMV+" = "#b256f5"), age=c("Young"="#2dcc30", "Old"="#f5ef38")),
  gap = unit(c(1, 4), "mm"),
  annotation_name_gp= gpar(fontsize = 12)
  # annotation_name_rot = 90
)

p <- Heatmap(
  gene_counts,
  name = "gene counts z-score",
  right_annotation = row_ha,
  # row_order = rownames(gene_counts),
  # column_order = colnames(gene_counts),
  heatmap_legend_param = list(
    legend_height = unit(4, "cm"),
    title_position = "lefttop-rot",
    title_gp = gpar(fontsize = 11)#, fontface = "bold")
  ),
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  # column_title = "GO BP 'Response to virus'",
  # column_title_gp = gpar(fontface = "bold", fontsize = 18)
)
p
setwd("/n/groups/kirschner/Will/HCMV/AdiposeTissue/deseq2")
png(file="heatmap.png", res = 400, width = 12, height = 4, units = 'in')
draw(p, column_title="GO BP 'Response to virus'",
     column_title_gp=grid::gpar(fontsize=16))
dev.off()

################################################################################
################################################################################
##### Virus
common_genes <- common_genes[common_genes %in% rownames(counts)]
gene_counts <- counts[common_genes, samples$participant_id]
gene_counts <- t(as.matrix(gene_counts))
gene_counts <- scale(gene_counts)
colnames(gene_counts) <- genes_unique[colnames(gene_counts), "entrezgene_accession"]
pathway <- "KEGG virus infections"
Heatmap(
  gene_counts,
  name = pathway,
  right_annotation = row_ha,
  # row_order = rownames(gene_counts),
  # column_order = colnames(gene_counts),
  heatmap_legend_param = list(
    legend_height = unit(4, "cm"),
    title_position = "lefttop-rot",
    title_gp = gpar(fontsize = 6, fontface = "bold")
  ),
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8)
)
################################################################################
### GO MF
library(org.Hs.eg.db)

ego <- gseGO(geneList     = deseq2_ranks,
             OrgDb        = org.Hs.eg.db,
             ont          = "MF",
             minGSSize    = 15,
             maxGSSize    = 1500,
             pvalueCutoff = 1,
             verbose      = TRUE,
             eps=0
)
ego_age <- gseGO(geneList     = deseq2_ranks_age,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "MF",
                 minGSSize    = 15,
                 maxGSSize    = 1500,
                 pvalueCutoff = 1,
                 verbose      = TRUE,
                 eps=0
)

ego_merged <- full_join(ego@result, ego_age@result, by = "ID", suffix = c("_virus","_age"))
ego_merged <- ego_merged[((!is.na(ego_merged$p.adjust_virus) & (ego_merged$p.adjust_virus < 0.05))) | (!is.na(ego_merged$p.adjust_age) & (ego_merged$p.adjust_age < 0.05)),]
ego_merged$Description <- ego_merged$Description_age
ego_merged$Description[is.na(ego_merged$Description)] <- ego_merged$Description_virus[is.na(ego_merged$Description)]
ego_merged$setSize <- ego_merged$setSize_age
ego_merged$setSize[is.na(ego_merged$setSize)] <- ego_merged$setSize_virus[is.na(ego_merged$setSize)]

ego_merged$core_enrichment_virus <- strsplit(ego_merged$core_enrichment_virus, "/")
ego_merged$core_enrichment_age <- strsplit(ego_merged$core_enrichment_age, "/")

pathway <- "structural constituent of ribosome"
n_top = 125

col <- "core_enrichment_age"

pathway_genes <- ego_merged[ego_merged$Description == pathway,][col][[1]][[1]][1:n_top]
pathway_genes <- pathway_genes[pathway_genes %in% rownames(counts)]
gene_counts <- counts[pathway_genes, samples$participant_id]
gene_counts <- t(as.matrix(gene_counts))
gene_counts <- scale(gene_counts)
colnames(gene_counts) <- genes_unique[colnames(gene_counts), "hgnc_symbol"]

Heatmap(
  gene_counts,
  name = pathway,
  right_annotation = row_ha,
  # row_order = rownames(gene_counts),
  # column_order = colnames(gene_counts),
  heatmap_legend_param = list(
    legend_height = unit(4, "cm"),
    title_position = "lefttop-rot",
    title_gp = gpar(fontsize = 6, fontface = "bold")
  ),
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8)
)

################################################################################
### GO CC
library(org.Hs.eg.db)

ego <- gseGO(geneList     = deseq2_ranks,
             OrgDb        = org.Hs.eg.db,
             ont          = "CC",
             minGSSize    = 15,
             maxGSSize    = 1500,
             pvalueCutoff = 1,
             verbose      = TRUE,
             eps=0
)
ego_age <- gseGO(geneList     = deseq2_ranks_age,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "CC",
                 minGSSize    = 15,
                 maxGSSize    = 1500,
                 pvalueCutoff = 1,
                 verbose      = TRUE,
                 eps=0
)

ego_merged <- full_join(ego@result, ego_age@result, by = "ID", suffix = c("_virus","_age"))
ego_merged <- ego_merged[((!is.na(ego_merged$p.adjust_virus) & (ego_merged$p.adjust_virus < 0.05))) | (!is.na(ego_merged$p.adjust_age) & (ego_merged$p.adjust_age < 0.05)),]
ego_merged$Description <- ego_merged$Description_age
ego_merged$Description[is.na(ego_merged$Description)] <- ego_merged$Description_virus[is.na(ego_merged$Description)]
ego_merged$setSize <- ego_merged$setSize_age
ego_merged$setSize[is.na(ego_merged$setSize)] <- ego_merged$setSize_virus[is.na(ego_merged$setSize)]

ego_merged$core_enrichment_virus <- strsplit(ego_merged$core_enrichment_virus, "/")
ego_merged$core_enrichment_age <- strsplit(ego_merged$core_enrichment_age, "/")

pathway <- "collagen-containing extracellular matrix"
n_top = 125

col <- "core_enrichment_virus"

pathway_genes <- ego_merged[ego_merged$Description == pathway,][col][[1]][[1]][1:n_top]
pathway_genes <- pathway_genes[pathway_genes %in% rownames(counts)]
gene_counts <- counts[pathway_genes, samples$participant_id]
gene_counts <- t(as.matrix(gene_counts))
gene_counts <- scale(gene_counts)
colnames(gene_counts) <- genes_unique[colnames(gene_counts), "hgnc_symbol"]

Heatmap(
  gene_counts,
  name = pathway,
  right_annotation = row_ha,
  # row_order = rownames(gene_counts),
  # column_order = colnames(gene_counts),
  heatmap_legend_param = list(
    legend_height = unit(4, "cm"),
    title_position = "lefttop-rot",
    title_gp = gpar(fontsize = 6, fontface = "bold")
  ),
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8)
)


################################################################################
### KEGG

kk <- gseKEGG(geneList     = deseq2_ranks,
              organism     = 'hsa',
              minGSSize    = 15,
              maxGSSize    = 1500,
              pvalueCutoff = 1,
              verbose      = TRUE,
              eps=0
)

kk_age <- gseKEGG(geneList     = deseq2_ranks_age,
                  organism     = 'hsa',
                  minGSSize    = 15,
                  maxGSSize    = 1500,
                  pvalueCutoff = 1,
                  verbose      = TRUE,
                  eps=0
)

kk_merged <- full_join(kk@result, kk_age@result, by = "ID", suffix = c("_virus","_age"))
kk_merged <- kk_merged[((!is.na(kk_merged$p.adjust_virus) & (kk_merged$p.adjust_virus < 0.05))) | (!is.na(kk_merged$p.adjust_age) & (kk_merged$p.adjust_age < 0.05)),]
kk_merged$Description <- kk_merged$Description_age
kk_merged$Description[is.na(kk_merged$Description)] <- kk_merged$Description_virus[is.na(kk_merged$Description)]
kk_merged$setSize <- kk_merged$setSize_age
kk_merged$setSize[is.na(kk_merged$setSize)] <- kk_merged$setSize_virus[is.na(kk_merged$setSize)]

kk_merged$core_enrichment_virus <- strsplit(kk_merged$core_enrichment_virus, "/")
kk_merged$core_enrichment_age <- strsplit(kk_merged$core_enrichment_age, "/")

pathway <- "Citrate cycle (TCA cycle)"
short_name <- pathway
n_top = 125

col <- "core_enrichment_virus"

pathway_genes <- kk_merged[kk_merged$Description == pathway,][col][[1]][[1]][1:n_top]
pathway_genes <- pathway_genes[pathway_genes %in% rownames(counts)]
gene_counts <- counts[pathway_genes, samples$participant_id]
gene_counts <- t(as.matrix(gene_counts))
gene_counts <- scale(gene_counts)
colnames(gene_counts) <- genes_unique[colnames(gene_counts), "hgnc_symbol"]

Heatmap(
  gene_counts,
  name = short_name,
  right_annotation = row_ha,
  # row_order = rownames(gene_counts),
  # column_order = colnames(gene_counts),
  heatmap_legend_param = list(
    legend_height = unit(4, "cm"),
    title_position = "lefttop-rot",
    title_gp = gpar(fontsize = 7, fontface = "bold")
  ),
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8)
)
