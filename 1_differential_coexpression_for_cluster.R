
# Differential co-expression analysis (DGCA + permutation)
#
# For each TCA gene × other gene pair, tests whether the Spearman correlation
# changes between Pre and Post bed rest using DGCA's permutation-based test.
#
# Why permutation instead of Fisher's z:
#   Fisher's z assumes N(0,1) for the test statistic, derived for Pearson r.
#   For Spearman r at n≈20, both the variance formula and the normality are
#   poor approximations. Permutation makes no distributional assumption:
#   it shuffles condition labels, recomputes z_diff, and the p-value is just
#   the empirical tail probability.
#
# Outputs:
#   dgca_all_tca_pairs.csv   — all TCA × other_gene results combined
#   <ENSEMBL>_dgca.csv       — one file per TCA gene, sorted by adj. p-value
#   tca_vs_insulin_perm.csv  — TCA gene × insulin permutation test (phenotype,
#                              handled separately since DGCA is gene-gene only)


if (!requireNamespace("DGCA", quietly = TRUE))
  stop("Install DGCA: install.packages('DGCA')")

suppressPackageStartupMessages({
  library(DGCA)
  library(data.table)
})

message("Loading counts (TPM)...")

# Helper: SRR_GSM_A1_POST_Homo_sapiens_RNA-Seq  →  A1_Post
tpm_short_name <- function(x) {
  p <- strsplit(x, "_")[[1]]
  paste0(p[3], "_",
         paste0(toupper(substr(p[4], 1, 1)), tolower(substr(p[4], 2, nchar(p[4])))))
}

tpm_raw_fp <- '/n/groups/kirschner/Will/BedRest/counts/rsem.merged.gene_tpm.tsv'
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


# ---- Parameters ----
N_HVG  <- 10000   # other-gene candidates: top N most variable genes
               # runtime scales as N_HVG² × N_PERM — this is the main lever
N_PERM <- 500    # permutations (use 1000 for final analysis)
OUT_DIR <- "/n/groups/kirschner/Will/BedRest/outputs/1_differential_coexpression"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

r_to_z <- function(r) 0.5 * log((1 + r) / (1 - r))
clip   <- function(r) pmax(pmin(r, 0.9999), -0.9999)

# =====================================================================
# PART 1: Gene–gene differential co-expression via DGCA
# =====================================================================

# DGCA design matrix: rows = samples, one binary column per condition
design_mat <- matrix(
  c(as.integer(grepl("_Post$", colnames(counts))),
    as.integer(grepl("_Pre$",  colnames(counts)))),
  ncol = 2,
  dimnames = list(colnames(counts), c("Post", "Pre"))
)

# Reduce to TCA genes + top variable other genes (runtime scales with n_genes²)
gene_vars <- apply(counts[other_genes, ], 1, var)
top_hvg   <- names(sort(gene_vars, decreasing = TRUE))[1:min(N_HVG, length(gene_vars))]
input_mat <- as.matrix(counts[c(tca_genes, top_hvg), ])

message(sprintf(
  "DGCA: %d genes × %d samples, %d permutations...",
  nrow(input_mat), ncol(input_mat), N_PERM
))

set.seed(42)
dgca_res <- ddcorAll(
  inputMat = input_mat,
  design   = design_mat,
  compare  = c("Post", "Pre"),   # cor1 = Post, cor2 = Pre
  corrType = "pearson",   # log2(TPM+1) is approx. normal so Pearson is appropriate.
                          # Also substantially faster than Spearman: Pearson uses
                          # BLAS-optimised matrix multiply; Spearman's ranking step
                          # breaks out of the BLAS path entirely.
  nPerm    = N_PERM,
  classify = TRUE,
  adjust   = "BH"
)


fwrite(as.data.table(dgca_res), file.path(OUT_DIR, "dgca_all.csv"))

# Column holding the (permutation-based) p-value differs by DGCA version
pval_col <- intersect(c("empPVal", "pValDiff"), colnames(dgca_res))[1]
padj_col <- intersect(c("pValDiff_adj", "empPVal_adj"), colnames(dgca_res))[1]
message(sprintf("Using p-value column: '%s', adjusted: '%s'", pval_col, padj_col))

# Keep only TCA × other_gene pairs; ensure TCA gene is always in Gene1
tca_g1 <- dgca_res[dgca_res$Gene1 %in% tca_genes & !(dgca_res$Gene2 %in% tca_genes), ]

tca_g2_raw <- dgca_res[dgca_res$Gene2 %in% tca_genes & !(dgca_res$Gene1 %in% tca_genes), ]
tca_g2 <- tca_g2_raw
tca_g2$Gene1 <- tca_g2_raw$Gene2   # cor1/cor2 are symmetric, no change needed
tca_g2$Gene2 <- tca_g2_raw$Gene1

tca_pairs <- rbind(tca_g1, tca_g2)
tca_pairs <- tca_pairs[order(tca_pairs[[padj_col]]), ]

fwrite(as.data.table(tca_pairs), file.path(OUT_DIR, "dgca_all_tca_pairs.csv"))

for (tca in tca_genes) {
  sub <- tca_pairs[tca_pairs$Gene1 == tca, ]
  fwrite(as.data.table(sub[order(sub[[padj_col]]), ]),
         file.path(OUT_DIR, paste0(tca, "_dgca.csv")))
}

message(sprintf("DGCA done. Significant pairs (adj.p<0.1): %d",
                sum(tca_pairs[[padj_col]] < 0.1, na.rm = TRUE)))


message(sprintf("Done. Output: %s", OUT_DIR))
