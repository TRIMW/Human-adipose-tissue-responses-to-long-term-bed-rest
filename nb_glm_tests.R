library(biomaRt)
library(msigdb)
library(dplyr)
library(ggplot2)


raw_counts_fp <- '/Users/willtrim/Documents/projs/bedrest/outputs/rsem.merged.gene_counts.tsv'
raw_counts <- read.csv(raw_counts_fp, sep="\t",row.names = 1)

raw_counts <- raw_counts[,!grepl("I2_POST", colnames(raw_counts))]

raw_counts$transcript_id.s. <- NULL

nrow(raw_counts)
min_counts <- 10
raw_counts <- raw_counts[rowSums(raw_counts) > min_counts,]
nrow(raw_counts)


counts_fp <- '/Users/willtrim/Documents/projs/bedrest/outputs/rlog_counts_star_rsem.csv'
counts <- read.csv(counts_fp, sep=",", row.names = 1)
counts <- counts[,!grepl("I2_Post", colnames(raw_counts))]
counts <- counts[rownames(counts) %in% rownames(raw_counts),]

gene_ids <- sapply(strsplit(rownames(counts), "\\."), head, 1)
rownames(counts) <- gene_ids

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

# counts$ensembl_gene_id <- gene_ids
# nrow(counts)
# counts <- merge(genes_unique, counts, by = "ensembl_gene_id")
# nrow(counts)
# rownames(counts) <- counts$entrezgene_id
# counts$ensembl_gene_id <- NULL
# counts$entrezgene_id <- NULL
# counts$X <- NULL
# 
# 
# colnames(counts) <- lapply(strsplit(colnames(counts), split = "_"), function(x) paste(x[[3]], x[[4]], sep="_"))

########
# phenotypes
pheno_fp <- '/Users/willtrim/Documents/projs/bedrest/outputs/BedRest_blood_phenotype.csv'
pheno <- read.csv(pheno_fp, sep=",")#, row.names = 1)


########
# KEGG
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

########
# I. Post Arterial Insulin vs TCA genes
########

# TCA genes post
pathway <- "Oxidative phosphorylation"
int_genes <- kegg_pathways[[pathway]]
int_genes_ensembl <- genes_unique[genes_unique$entrezgene_id %in% int_genes, ] %>% pull(ensembl_gene_id)

int_counts <- t(counts[rownames(counts) %in% int_genes_ensembl, grepl("Post", colnames(counts))])
other_counts <- t(counts[!(rownames(counts) %in% int_genes_ensembl), grepl("Post", colnames(counts))])
post_counts <- t(counts[, grepl("Post", colnames(counts))])

rownames(int_counts) <- lapply(strsplit(rownames(int_counts), split = "_"), function(x) x[[1]])
rownames(other_counts) <- lapply(strsplit(rownames(other_counts), split = "_"), function(x) x[[1]])
rownames(post_counts) <- lapply(strsplit(rownames(post_counts), split = "_"), function(x) x[[1]])

# Arterial Insulin post
pheno_name <- "Art..Insulin..mIU.L."
post_pheno <- pheno[pheno$Time.Point == "Post Bed Rest", ]
int_pheno <- as.data.frame(post_pheno[,pheno_name])
rownames(int_pheno) <- post_pheno$Participant.ID
colnames(int_pheno) <- c(pheno_name)


predictors <- int_genes_ensembl[int_genes_ensembl %in% rownames(counts)]#[1]
other_genes <- rownames(counts)[!(rownames(counts) %in% int_genes_ensembl)]
# 
# formula_obj <- reformulate(predictors, response = pheno_name)



int_gene_solo_results <- list()


counts_pheno <- merge(post_counts, int_pheno, by = "row.names")
rownames(counts_pheno) <- counts_pheno$Row.names
counts_pheno$Row.names <- NULL


for(int_gene in predictors) {
  
  formula_obj <- reformulate(c(int_gene), response = pheno_name)
  
  base_model <- glm(formula_obj,
                    family = Gamma(link = "log"), data = counts_pheno)
  
  
  results <- data.frame(
    gene = character(),
    lrt_pvalue = numeric(),
    delta_aic = numeric()
  )
  
  for(other_gene in other_genes) {
    formula_obj <- reformulate(c(int_gene, other_gene), response = pheno_name)
    
    test_model <- glm(formula_obj,
                      family = Gamma(link = "log"), data = counts_pheno)
    
    lrt <- anova(base_model, test_model, test = "Chisq")
    
    results <- rbind(results, data.frame(
      gene = other_gene,
      lrt_pvalue = lrt$`Pr(>Chi)`[2],
      delta_aic = AIC(base_model) - AIC(test_model)
    ))
    # break
  }
  
  int_gene_solo_results[[int_gene]] <- results
}



setwd("/Users/willtrim/Documents/projs/bedrest/outputs/ox_phos_solo_insulin_art_post")
library(data.table)

for (int_gene in names(int_gene_solo_results)) {
  fwrite(
    int_gene_solo_results[[int_gene]],
    paste(int_gene, ".csv", sep=""),
    row.names = FALSE
  )
}

