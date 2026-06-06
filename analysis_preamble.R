
# Shared data loading for all analysis scripts.
# Source this once per R session: source("analysis_preamble.R")
# Each analysis script calls `if (!exists("counts")) source("analysis_preamble.R")`

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(purrr)
  library(data.table)
  library(clusterProfiler)
})

message("Loading counts (TPM)...")

# Helper: SRR_GSM_A1_POST_Homo_sapiens_RNA-Seq  →  A1_Post
tpm_short_name <- function(x) {
  p <- strsplit(x, "_")[[1]]
  paste0(p[3], "_",
         paste0(toupper(substr(p[4], 1, 1)), tolower(substr(p[4], 2, nchar(p[4])))))
}

tpm_raw_fp <- '/Users/willtrim/Documents/projs/bedrest/outputs/rsem.merged.gene_tpm.tsv'
tpm_raw    <- read.csv(tpm_raw_fp, sep="\t", row.names=1, check.names=FALSE)
tpm_raw[["transcript_id(s)"]] <- NULL
tpm_raw    <- tpm_raw[, !grepl("I2_POST", colnames(tpm_raw))]          # exclude outlier
colnames(tpm_raw) <- sapply(colnames(tpm_raw), tpm_short_name)          # A1_Post, A1_Pre …
rownames(tpm_raw) <- sapply(strsplit(rownames(tpm_raw), "\\."), `[[`, 1) # strip ENSEMBL version

# Filter: keep genes with mean raw TPM >= 1 across all retained samples.
# Genes below this threshold are dominated by Poisson counting noise;
# log-transform does not stabilise their variance.
keep    <- rowMeans(tpm_raw) >= 1
tpm_raw <- tpm_raw[keep, ]
message(sprintf("Genes after mean TPM >= 1 filter: %d", nrow(tpm_raw)))

# Log2(TPM + 1): stabilises NB overdispersion for the retained expressed genes.
counts <- log2(tpm_raw + 1)

# raw_counts: kept in memory only for script 6 (MuSiC deconvolution needs raw integers).
raw_counts_fp <- '/Users/willtrim/Documents/projs/bedrest/outputs/rsem.merged.gene_counts.tsv'
raw_counts    <- read.csv(raw_counts_fp, sep="\t", row.names=1, check.names=FALSE)
raw_counts[["transcript_id(s)"]] <- NULL
raw_counts    <- raw_counts[, !grepl("I2_POST", colnames(raw_counts))]
colnames(raw_counts) <- sapply(colnames(raw_counts), tpm_short_name)
rownames(raw_counts) <- sapply(strsplit(rownames(raw_counts), "\\."), `[[`, 1)
raw_counts    <- raw_counts[rownames(raw_counts) %in% rownames(counts), ]

pheno_fp <- '/Users/willtrim/Documents/projs/bedrest/outputs/BedRest_blood_phenotype.csv'
pheno    <- read.csv(pheno_fp)

# Phenotype column names after read.csv
pheno_vars <- c(
  "Venous.Glucose..mmol.L.", "Venous.Insulin..mU.L.",
  "Art..Glucose..mmol.L.",   "Fat.ox.kg",        "cho.ox.Kg",
  "FFA..mmol.L.",            "TAG..mmol.L.",
  "Art..Insulin..mIU.L.",    "Glucose.Disposal",
  "Post.clamp.Fat.ox.kg", "Post.clamp.cho.ox.Kg"
)

# ---- BiomaRt: ENSEMBL → Entrez (for KEGG) and ENSEMBL → Symbol ----
message("Querying BiomaRt (this can take ~1 min)...")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

genes_entrez <- getBM(
  filters    = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  values     = rownames(counts), mart = mart
)
genes_unique <- genes_entrez %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  group_by(entrezgene_id) %>% filter(n() == 1) %>% ungroup()

sym_map <- getBM(
  filters    = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values     = rownames(counts), mart = mart
)
sym_map <- sym_map %>%
  filter(hgnc_symbol != "") %>%
  group_by(hgnc_symbol) %>% filter(n() == 1) %>% ungroup() %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Gene-symbol keyed expression matrix (for tools that need symbols)
counts_sym           <- counts[sym_map$ensembl_gene_id, ]
rownames(counts_sym) <- sym_map$hgnc_symbol

# ---- KEGG pathways ----
message("Downloading KEGG pathways...")
kegg_df   <- download_KEGG(species = "hsa")
kegg_list <- kegg_df[["KEGGPATHID2EXTID"]] %>% group_by(from) %>% group_split()
names_map <- setNames(kegg_df[["KEGGPATHID2NAME"]]$to,
                      kegg_df[["KEGGPATHID2NAME"]]$from)
kegg_sets            <- map(kegg_list, ~pull(.x, to))
names(kegg_sets)     <- map_chr(kegg_list, ~names_map[.x[[1, 1]]])

# ---- Pathway selection ----
# Change this to "Oxidative phosphorylation" if needed.
# Run grep("citrate|TCA", names(kegg_sets), ignore.case=TRUE) to verify name.
pathway <- "Citrate cycle (TCA cycle)"

int_genes_entrez <- kegg_sets[[pathway]]
tca_ensembl      <- genes_unique %>%
  filter(entrezgene_id %in% int_genes_entrez) %>%
  pull(ensembl_gene_id)
tca_genes   <- tca_ensembl[tca_ensembl %in% rownames(counts)]
other_genes <- rownames(counts)[!(rownames(counts) %in% tca_genes)]

message(sprintf("Pathway '%s': %d genes in data", pathway, length(tca_genes)))

# ---- Split matrices ----
pre_mat  <- t(counts[, grep("_Pre$",  colnames(counts), value = TRUE)])
post_mat <- t(counts[, grep("_Post$", colnames(counts), value = TRUE)])
rownames(pre_mat)  <- sub("_Pre$",  "", rownames(pre_mat))
rownames(post_mat) <- sub("_Post$", "", rownames(post_mat))

# Paired delta (only participants with both timepoints)
pre_ids  <- sub("_Pre$",  "", grep("_Pre$",  colnames(counts), value = TRUE))
post_ids <- sub("_Post$", "", grep("_Post$", colnames(counts), value = TRUE))
paired   <- intersect(pre_ids, post_ids)

pre_expr  <- counts[, paste0(paired, "_Pre")]
post_expr <- counts[, paste0(paired, "_Post")]
colnames(pre_expr) <- colnames(post_expr) <- paired
delta_expr <- post_expr - pre_expr   # genes × participants

# Delta phenotypes
pheno_post <- pheno[pheno$Time.Point == "Post Bed Rest", ]
pheno_pre  <- pheno[pheno$Time.Point == "Pre Bed Rest",  ]
rownames(pheno_post) <- pheno_post$Participant.ID
rownames(pheno_pre)  <- pheno_pre$Participant.ID
pheno_paired <- intersect(paired, rownames(pheno_post))

delta_pheno <- as.data.frame(
  pheno_post[pheno_paired, pheno_vars] - pheno_pre[pheno_paired, pheno_vars]
)

message("Preamble complete.")
