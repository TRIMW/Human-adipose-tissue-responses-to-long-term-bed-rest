
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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

if (!requireNamespace("DGCA", quietly = TRUE))
  stop("Install DGCA: install.packages('DGCA')")

suppressPackageStartupMessages({
  library(DGCA)
  library(data.table)
})
if (!exists("counts")) source("analysis_preamble.R")

# ---- Parameters ----
N_HVG  <- 5000   # other-gene candidates: top N most variable genes
               # runtime scales as N_HVG² × N_PERM — this is the main lever
N_PERM <- 250    # permutations (use 1000 for final analysis)
OUT_DIR <- "/Users/willtrim/Documents/projs/bedrest/outputs/1_differential_coexpression"
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

# =====================================================================
# PART 2: TCA gene × arterial insulin  (manual permutation)
# DGCA only handles gene–gene pairs; phenotype correlations done here.
# =====================================================================

# Attach the correct insulin value to every sample
all_cols  <- colnames(counts)
is_post   <- grepl("_Post$", all_cols)
part_ids  <- sub("_(Pre|Post)$", "", all_cols)

rownames(pheno_post) <- pheno_post$Participant.ID
rownames(pheno_pre)  <- pheno_pre$Participant.ID

insulin_all <- ifelse(
  is_post,
  pheno_post[part_ids, "Art..Insulin..mIU.L."],
  pheno_pre[ part_ids, "Art..Insulin..mIU.L."]
)
names(insulin_all) <- all_cols

ok <- !is.na(insulin_all)   # samples with insulin data

perm_one_tca <- function(expr_row) {
  # Observed: correlation in Post vs Pre with insulin
  r_post <- cor(expr_row[is_post & ok],  insulin_all[is_post & ok],
                method = "pearson", use = "complete.obs")
  r_pre  <- cor(expr_row[!is_post & ok], insulin_all[!is_post & ok],
                method = "pearson", use = "complete.obs")
  obs    <- r_to_z(clip(r_post)) - r_to_z(clip(r_pre))

  # Permute condition labels (reshuffle which samples are "post")
  perm_diffs <- replicate(N_PERM, {
    perm_post <- sample(is_post)
    rp <- cor(expr_row[perm_post & ok],  insulin_all[perm_post & ok],
              method = "pearson", use = "complete.obs")
    rr <- cor(expr_row[!perm_post & ok], insulin_all[!perm_post & ok],
              method = "pearson", use = "complete.obs")
    r_to_z(clip(rp)) - r_to_z(clip(rr))
  })

  p_val <- mean(abs(perm_diffs) >= abs(obs))

  data.frame(r_pre = r_pre, r_post = r_post,
             z_diff = obs, perm_p = p_val)
}

message("Permutation test: TCA genes × arterial insulin...")
insulin_res <- do.call(rbind, lapply(tca_genes, function(tca) {
  res <- perm_one_tca(as.numeric(counts[tca, ]))
  cbind(tca_gene = tca, res)
}))
insulin_res$perm_p_adj <- p.adjust(insulin_res$perm_p, method = "BH")
insulin_res <- insulin_res[order(insulin_res$perm_p), ]

fwrite(as.data.table(insulin_res),
       file.path(OUT_DIR, "tca_vs_insulin_perm.csv"))

message(sprintf("Done. Output: %s", OUT_DIR))
